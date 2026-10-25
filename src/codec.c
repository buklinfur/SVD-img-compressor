#include "codec.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// .svdc layout, version 1 (all integers little-endian, floats IEEE-754 binary32):
//
//   offset  size  field
//   0       4     magic "SVDC"
//   4       1     version (1)
//   5       1     bits per factor value: 8, 16 or 32
//   6       1     channels (1..4)
//   7       1     reserved (0)
//   8       4     width  (u32)
//   12      4     height (u32)
//   then, for every channel:
//     u8  type: 0 = constant, 1 = low rank
//     constant: u8 value
//     low rank: u32 k, f32 s[k],
//               k columns of U: [f32 scale] + height values,
//               k columns of V: [f32 scale] + width values
//   a value is int8 (q / 127 * scale), int16 (q / 32767 * scale) or a raw f32;
//   the per-column scale (max |entry|) is present for 8 and 16 bits only

#define SVDC_VERSION 1
#define SVDC_HEADER_SIZE 16
#define SVDC_MAX_SIDE (1 << 24)
#define SVDC_MAX_PIXELS ((uint64_t) 1 << 30)

_Static_assert (sizeof (float) == 4, "the .svdc format stores IEEE-754 binary32 floats");

CompressOptions defaultCompressOptions (void) {
    CompressOptions opt;
    opt.k = 50;
    opt.energy = 0.0;
    opt.bits = 8;
    opt.threads = 0;
    opt.svd = defaultSvdOptions ();
    return opt;
}

const char* codecErrorString (int code) {
    switch (code) {
        case SVDC_OK: return "ok";
        case SVDC_ERR_ALLOC: return "out of memory";
        case SVDC_ERR_IO: return "file read/write error";
        case SVDC_ERR_FORMAT: return "not an .svdc file or unsupported version";
        case SVDC_ERR_TRUNCATED: return "file is truncated";
        case SVDC_ERR_CORRUPT: return "file is corrupt";
        case SVDC_ERR_ARGUMENT: return "invalid argument";
        default: return "unknown error";
    }
}

void freeCompressed (Compressed* c) {
    if (c == NULL)
        return;
    for (int ch = 0; ch < 4; ch++)
        freeSvd (c->channel[ch].svd);
    free (c);
}

// ---------------------------------------------------------------------------
// compression

static int compressChannel (const Image* img, int ch, const CompressOptions* opt, CompressedChannel* out) {
    unsigned char value = 0;
    if (isConstantChannel (img, ch, &value)) {
        out->constant = 1;
        out->value = value;
        return SVDC_OK;
    }

    Matrix* a = channelToMatrix (img, ch);
    if (a == NULL)
        return SVDC_ERR_ALLOC;
    SvdOptions svdOpt = opt->svd;
    svdOpt.seed += (uint64_t) ch;  // independent sketches, same result for any thread count
    out->svd = opt->energy > 0.0 ? svdForEnergy (a, opt->energy, &svdOpt)
                                 : svdCompute (a, opt->k, &svdOpt);
    freeMatrix (a);
    return out->svd != NULL ? SVDC_OK : SVDC_ERR_ALLOC;
}

Compressed* compressImage (const Image* img, const CompressOptions* opt) {
    const CompressOptions defaults = defaultCompressOptions ();
    if (opt == NULL)
        opt = &defaults;
    if (opt->bits != 8 && opt->bits != 16 && opt->bits != 32)
        return NULL;

    Compressed* c = (Compressed*) calloc (1, sizeof (Compressed));
    if (c == NULL)
        return NULL;
    c->width = img->width;
    c->height = img->height;
    c->channels = img->channels;
    c->bits = opt->bits;

    int failed = 0;
#ifdef _OPENMP
    int threads = opt->threads > 0 ? opt->threads : omp_get_max_threads ();
    if (threads > img->channels)
        threads = img->channels;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(| : failed)
#endif
    for (int ch = 0; ch < img->channels; ch++)
        failed |= compressChannel (img, ch, opt, &c->channel[ch]) != SVDC_OK;

    if (failed) {
        freeCompressed (c);
        return NULL;
    }
    return c;
}

Image* decompressImage (const Compressed* c) {
    Image* img = makeImage (c->width, c->height, c->channels);
    if (img == NULL)
        return NULL;
    const size_t pixels = (size_t) c->width * (size_t) c->height;
    const size_t stride = (size_t) c->channels;

    for (int ch = 0; ch < c->channels; ch++) {
        const CompressedChannel* cc = &c->channel[ch];
        if (cc->constant || cc->svd == NULL) {
            for (size_t i = 0; i < pixels; i++)
                img->data[i * stride + (size_t) ch] = cc->value;
            continue;
        }
        Matrix* m = svdReconstruct (cc->svd, 0);
        if (m == NULL) {
            freeImage (img);
            return NULL;
        }
        matrixToChannel (m, img, ch);
        freeMatrix (m);
    }
    return img;
}

// ---------------------------------------------------------------------------
// serialization helpers

static int quantLevels (int bits) {
    return bits == 8 ? 127 : 32767;
}

static size_t valueBytes (int bits) {
    return (size_t) bits / 8;
}

// bytes of one factor column holding `len` values
static size_t columnBytes (int bits, size_t len) {
    return (bits == 32 ? 0 : 4) + len * valueBytes (bits);
}

static size_t channelBytes (const CompressedChannel* cc, int bits, int width, int height) {
    if (cc->constant || cc->svd == NULL)
        return 2;
    const size_t k = (size_t) cc->svd->rank;
    return 1 + 4 + 4 * k + k * columnBytes (bits, (size_t) height) + k * columnBytes (bits, (size_t) width);
}

size_t encodedSize (const Compressed* c) {
    size_t size = SVDC_HEADER_SIZE;
    for (int ch = 0; ch < c->channels; ch++)
        size += channelBytes (&c->channel[ch], c->bits, c->width, c->height);
    return size;
}

typedef struct Writer {
    unsigned char* p;
    size_t pos;
} Writer;

static void putU8 (Writer* w, unsigned value) {
    w->p[w->pos++] = (unsigned char) (value & 0xFFu);
}

static void putU16 (Writer* w, unsigned value) {
    putU8 (w, value);
    putU8 (w, value >> 8);
}

static void putU32 (Writer* w, uint32_t value) {
    for (int i = 0; i < 4; i++)
        putU8 (w, (unsigned) (value >> (8 * i)));
}

static void putF32 (Writer* w, double value) {
    float f = (float) value;
    uint32_t bitsValue;
    memcpy (&bitsValue, &f, sizeof bitsValue);
    putU32 (w, bitsValue);
}

// write column r of m (all rows) with the requested precision
static void putColumn (Writer* w, const Matrix* m, int r, int bits) {
    if (bits == 32) {
        for (int i = 0; i < m->rows; i++)
            putF32 (w, matrixGet (m, i, r));
        return;
    }
    double scale = 0.0;
    for (int i = 0; i < m->rows; i++) {
        double v = fabs (matrixGet (m, i, r));
        if (v > scale)
            scale = v;
    }
    // store the scale as float and quantize against that exact float value,
    // so that the decoder sees the same scale the encoder used
    const double storedScale = (double) (float) scale;
    putF32 (w, storedScale);
    const int levels = quantLevels (bits);
    for (int i = 0; i < m->rows; i++) {
        long q = 0;
        if (storedScale > 0.0) {
            q = lround (matrixGet (m, i, r) / storedScale * levels);
            if (q > levels)
                q = levels;
            if (q < -levels)
                q = -levels;
        }
        // two's complement encoding of the signed value
        unsigned u = (unsigned) (q < 0 ? q + (bits == 8 ? 256L : 65536L) : q);
        if (bits == 8)
            putU8 (w, u);
        else
            putU16 (w, u);
    }
}

int encodeCompressed (const Compressed* c, unsigned char** out, size_t* size) {
    if (c->bits != 8 && c->bits != 16 && c->bits != 32)
        return SVDC_ERR_ARGUMENT;
    const size_t total = encodedSize (c);
    Writer w;
    w.p = (unsigned char*) malloc (total);
    w.pos = 0;
    if (w.p == NULL)
        return SVDC_ERR_ALLOC;

    putU8 (&w, 'S');
    putU8 (&w, 'V');
    putU8 (&w, 'D');
    putU8 (&w, 'C');
    putU8 (&w, SVDC_VERSION);
    putU8 (&w, (unsigned) c->bits);
    putU8 (&w, (unsigned) c->channels);
    putU8 (&w, 0);
    putU32 (&w, (uint32_t) c->width);
    putU32 (&w, (uint32_t) c->height);

    for (int ch = 0; ch < c->channels; ch++) {
        const CompressedChannel* cc = &c->channel[ch];
        if (cc->constant || cc->svd == NULL) {
            putU8 (&w, 0);
            putU8 (&w, cc->value);
            continue;
        }
        const Svd* svd = cc->svd;
        putU8 (&w, 1);
        putU32 (&w, (uint32_t) svd->rank);
        for (int r = 0; r < svd->rank; r++)
            putF32 (&w, svd->s[r]);
        for (int r = 0; r < svd->rank; r++)
            putColumn (&w, svd->u, r, c->bits);
        for (int r = 0; r < svd->rank; r++)
            putColumn (&w, svd->v, r, c->bits);
    }

    *out = w.p;
    *size = w.pos;
    return SVDC_OK;
}

// ---------------------------------------------------------------------------
// parsing

typedef struct Reader {
    const unsigned char* p;
    size_t size;
    size_t pos;
} Reader;

static int has (const Reader* r, size_t n) {
    return r->size - r->pos >= n;
}

static unsigned getU8 (Reader* r) {
    return r->p[r->pos++];
}

static unsigned getU16 (Reader* r) {
    unsigned lo = getU8 (r);
    return lo | (getU8 (r) << 8);
}

static uint32_t getU32 (Reader* r) {
    uint32_t v = 0;
    for (int i = 0; i < 4; i++)
        v |= (uint32_t) getU8 (r) << (8 * i);
    return v;
}

static double getF32 (Reader* r) {
    uint32_t bitsValue = getU32 (r);
    float f;
    memcpy (&f, &bitsValue, sizeof f);
    return (double) f;
}

// read one factor column into column `col` of m; the caller checked the size
static int getColumn (Reader* r, Matrix* m, int col, int bits) {
    if (bits == 32) {
        for (int i = 0; i < m->rows; i++) {
            double v = getF32 (r);
            if (!isfinite (v))
                return SVDC_ERR_CORRUPT;
            matrixSet (m, i, col, v);
        }
        return SVDC_OK;
    }
    const double scale = getF32 (r);
    if (!isfinite (scale) || scale < 0.0)
        return SVDC_ERR_CORRUPT;
    const double step = scale / quantLevels (bits);
    for (int i = 0; i < m->rows; i++) {
        long q = bits == 8 ? (long) getU8 (r) : (long) getU16 (r);
        if (bits == 8 && q > 127)
            q -= 256;
        else if (bits == 16 && q > 32767)
            q -= 65536;
        matrixSet (m, i, col, (double) q * step);
    }
    return SVDC_OK;
}

int decodeCompressed (const unsigned char* buf, size_t size, Compressed** out) {
    *out = NULL;
    Reader r = { buf, size, 0 };
    if (!has (&r, SVDC_HEADER_SIZE))
        return size >= 4 && memcmp (buf, "SVDC", 4) == 0 ? SVDC_ERR_TRUNCATED : SVDC_ERR_FORMAT;
    if (memcmp (buf, "SVDC", 4) != 0)
        return SVDC_ERR_FORMAT;
    r.pos = 4;
    const unsigned version = getU8 (&r);
    const unsigned bits = getU8 (&r);
    const unsigned channels = getU8 (&r);
    (void) getU8 (&r);  // reserved
    const uint32_t width = getU32 (&r);
    const uint32_t height = getU32 (&r);

    if (version != SVDC_VERSION)
        return SVDC_ERR_FORMAT;
    if ((bits != 8 && bits != 16 && bits != 32) || channels < 1 || channels > 4)
        return SVDC_ERR_CORRUPT;
    if (width == 0 || height == 0 || width > SVDC_MAX_SIDE || height > SVDC_MAX_SIDE ||
        (uint64_t) width * height > SVDC_MAX_PIXELS)
        return SVDC_ERR_CORRUPT;

    Compressed* c = (Compressed*) calloc (1, sizeof (Compressed));
    if (c == NULL)
        return SVDC_ERR_ALLOC;
    c->width = (int) width;
    c->height = (int) height;
    c->channels = (int) channels;
    c->bits = (int) bits;
    const uint32_t maxRank = width < height ? width : height;

    int status = SVDC_OK;
    for (int ch = 0; ch < c->channels && status == SVDC_OK; ch++) {
        CompressedChannel* cc = &c->channel[ch];
        if (!has (&r, 1)) {
            status = SVDC_ERR_TRUNCATED;
            break;
        }
        const unsigned type = getU8 (&r);
        if (type == 0) {
            if (!has (&r, 1)) {
                status = SVDC_ERR_TRUNCATED;
                break;
            }
            cc->constant = 1;
            cc->value = (unsigned char) getU8 (&r);
            continue;
        }
        if (type != 1) {
            status = SVDC_ERR_CORRUPT;
            break;
        }
        if (!has (&r, 4)) {
            status = SVDC_ERR_TRUNCATED;
            break;
        }
        const uint32_t k = getU32 (&r);
        if (k == 0 || k > maxRank) {
            status = SVDC_ERR_CORRUPT;
            break;
        }
        // check the payload size before allocating anything it describes
        const size_t need = 4 * (size_t) k + (size_t) k * columnBytes ((int) bits, height) +
                            (size_t) k * columnBytes ((int) bits, width);
        if (!has (&r, need)) {
            status = SVDC_ERR_TRUNCATED;
            break;
        }

        cc->svd = makeSvd ((int) height, (int) width, (int) k);
        if (cc->svd == NULL) {
            status = SVDC_ERR_ALLOC;
            break;
        }
        for (uint32_t i = 0; i < k; i++) {
            double s = getF32 (&r);
            if (!isfinite (s) || s < 0.0)
                status = SVDC_ERR_CORRUPT;
            cc->svd->s[i] = s;
        }
        for (int i = 0; i < (int) k && status == SVDC_OK; i++)
            status = getColumn (&r, cc->svd->u, i, (int) bits);
        for (int i = 0; i < (int) k && status == SVDC_OK; i++)
            status = getColumn (&r, cc->svd->v, i, (int) bits);
    }

    if (status == SVDC_OK && r.pos != r.size)
        status = SVDC_ERR_CORRUPT;  // trailing garbage
    if (status != SVDC_OK) {
        freeCompressed (c);
        return status;
    }
    *out = c;
    return SVDC_OK;
}

// ---------------------------------------------------------------------------
// files

int writeCompressedFile (const Compressed* c, const char* filename, size_t* written) {
    unsigned char* buf = NULL;
    size_t size = 0;
    int status = encodeCompressed (c, &buf, &size);
    if (status != SVDC_OK)
        return status;
    FILE* f = fopen (filename, "wb");
    if (f == NULL) {
        free (buf);
        return SVDC_ERR_IO;
    }
    size_t done = fwrite (buf, 1, size, f);
    int closed = fclose (f);
    free (buf);
    if (done != size || closed != 0)
        return SVDC_ERR_IO;
    if (written != NULL)
        *written = size;
    return SVDC_OK;
}

int readCompressedFile (const char* filename, Compressed** out, size_t* size) {
    *out = NULL;
    FILE* f = fopen (filename, "rb");
    if (f == NULL)
        return SVDC_ERR_IO;

    // read in chunks: works for pipes too and needs no ftell/fseek size tricks
    size_t capacity = 1 << 16, length = 0;
    unsigned char* buf = (unsigned char*) malloc (capacity);
    int status = buf != NULL ? SVDC_OK : SVDC_ERR_ALLOC;
    while (status == SVDC_OK) {
        if (length == capacity) {
            if (capacity > SIZE_MAX / 2) {
                status = SVDC_ERR_ALLOC;
                break;
            }
            unsigned char* bigger = (unsigned char*) realloc (buf, capacity * 2);
            if (bigger == NULL) {
                status = SVDC_ERR_ALLOC;
                break;
            }
            buf = bigger;
            capacity *= 2;
        }
        size_t got = fread (buf + length, 1, capacity - length, f);
        length += got;
        if (got == 0) {
            if (ferror (f))
                status = SVDC_ERR_IO;
            break;
        }
    }
    fclose (f);

    if (status == SVDC_OK)
        status = decodeCompressed (buf, length, out);
    if (status == SVDC_OK && size != NULL)
        *size = length;
    free (buf);
    return status;
}
