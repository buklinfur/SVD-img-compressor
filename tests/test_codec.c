// tests for image I/O, the .svdc format and the full compress/decompress pipeline

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "codec.h"
#include "image.h"
#include "metrics.h"
#include "test_util.h"

#ifndef EXAMPLES_DIR
#define EXAMPLES_DIR "examples"
#endif

// deterministic smooth test pattern with some texture
static Image* patternImage (int width, int height, int channels) {
    Image* img = makeImage (width, height, channels);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            for (int c = 0; c < channels; c++) {
                double v = 128.0 + 60.0 * sin (0.05 * x * (c + 1)) * cos (0.07 * y) +
                           30.0 * sin (0.3 * (x + y + 7 * c));
                img->data[((size_t) y * (size_t) width + (size_t) x) * (size_t) channels + (size_t) c] =
                    (unsigned char) (v < 0 ? 0 : v > 255 ? 255 : v);
            }
    return img;
}

static Image* roundTrip (const Image* img, const CompressOptions* opt, size_t* size) {
    Compressed* c = compressImage (img, opt);
    CHECK (c != NULL);
    unsigned char* buf = NULL;
    size_t len = 0;
    CHECK (encodeCompressed (c, &buf, &len) == SVDC_OK);
    CHECK (len == encodedSize (c));
    Compressed* d = NULL;
    CHECK (decodeCompressed (buf, len, &d) == SVDC_OK);
    Image* out = d != NULL ? decompressImage (d) : NULL;
    if (size != NULL)
        *size = len;
    free (buf);
    freeCompressed (c);
    freeCompressed (d);
    return out;
}

static void testMatrixToChannelClamps (void) {
    Image* img = makeImage (4, 1, 1);
    Matrix* m = makeMatrix (1, 4);
    m->data[0] = -3.0;
    m->data[1] = 254.6;
    m->data[2] = 1000.0;
    m->data[3] = (double) NAN;
    matrixToChannel (m, img, 0);
    CHECK (img->data[0] == 0);
    CHECK (img->data[1] == 255);
    CHECK (img->data[2] == 255);
    CHECK (img->data[3] == 0);
    freeMatrix (m);
    freeImage (img);
}

static void testFullRankIsLossless (void) {
    // exact SVD + float32 factors + all triplets: rounding restores every pixel
    Image* img = patternImage (37, 23, 3);
    CompressOptions opt = defaultCompressOptions ();
    opt.svd.method = SVD_JACOBI;
    opt.k = 23;
    opt.bits = 32;
    Image* out = roundTrip (img, &opt, NULL);
    CHECK (out != NULL && imagePsnr (img, out) == (double) INFINITY);
    freeImage (out);
    freeImage (img);
}

static void testQualityGrowsWithK (void) {
    Image* img = patternImage (120, 80, 3);
    CompressOptions opt = defaultCompressOptions ();
    double previous = 0.0;
    size_t previousSize = 0;
    const int ks[] = { 1, 4, 16, 40 };
    for (int i = 0; i < 4; i++) {
        opt.k = ks[i];
        size_t size = 0;
        Image* out = roundTrip (img, &opt, &size);
        double psnr = imagePsnr (img, out);
        CHECK (psnr > previous);
        CHECK (size > previousSize);
        previous = psnr;
        previousSize = size;
        freeImage (out);
    }
    CHECK (previous > 35.0);
    freeImage (img);
}

static void testQuantizationLoss (void) {
    // 16-bit factors should be indistinguishable from float32. 8-bit factors add
    // quantization noise that caps quality at roughly 50 dB, invisible in practice
    Image* img = patternImage (90, 70, 1);
    CompressOptions opt = defaultCompressOptions ();
    opt.k = 30;
    opt.bits = 32;
    Image* f32 = roundTrip (img, &opt, NULL);
    opt.bits = 16;
    Image* i16 = roundTrip (img, &opt, NULL);
    opt.bits = 8;
    Image* i8 = roundTrip (img, &opt, NULL);
    double p32 = imagePsnr (img, f32), p16 = imagePsnr (img, i16), p8 = imagePsnr (img, i8);
    CHECK (fabs (p32 - p16) < 0.05);
    CHECK (p8 > 45.0 && p8 <= p32 + 0.1);
    freeImage (f32);
    freeImage (i16);
    freeImage (i8);
    freeImage (img);
}

static void testConstantAndGrayChannels (void) {
    // gray + alpha where alpha is constant: stored as one byte, restored exactly
    Image* img = patternImage (50, 30, 2);
    for (int i = 0; i < 50 * 30; i++)
        img->data[2 * i + 1] = 200;
    CompressOptions opt = defaultCompressOptions ();
    opt.k = 10;
    Compressed* c = compressImage (img, &opt);
    CHECK (!c->channel[0].constant && c->channel[0].svd != NULL);
    CHECK (c->channel[1].constant && c->channel[1].value == 200 && c->channel[1].svd == NULL);
    Image* out = decompressImage (c);
    int alphaOk = 1;
    for (int i = 0; i < 50 * 30; i++)
        alphaOk &= out->data[2 * i + 1] == 200;
    CHECK (alphaOk);
    freeImage (out);
    freeCompressed (c);

    // a completely black image is all constant channels
    Image* black = makeImage (8, 8, 3);
    Image* restored = roundTrip (black, &opt, NULL);
    CHECK (imagePsnr (black, restored) == (double) INFINITY);
    freeImage (restored);
    freeImage (black);
    freeImage (img);
}

static void testEnergyMode (void) {
    Image* img = patternImage (100, 60, 1);
    CompressOptions opt = defaultCompressOptions ();
    opt.energy = 0.9999;
    Compressed* c = compressImage (img, &opt);
    CHECK (c != NULL && c->channel[0].svd != NULL);
    int k1 = c->channel[0].svd->rank;
    freeCompressed (c);
    opt.energy = 0.99;
    c = compressImage (img, &opt);
    int k2 = c->channel[0].svd->rank;
    CHECK (k2 >= 1 && k2 <= k1);
    freeCompressed (c);
    freeImage (img);
}

static void testRejectsBadInput (void) {
    Image* img = patternImage (16, 12, 3);
    CompressOptions opt = defaultCompressOptions ();
    opt.k = 4;
    Compressed* c = compressImage (img, &opt);
    unsigned char* buf = NULL;
    size_t len = 0;
    encodeCompressed (c, &buf, &len);
    Compressed* d = NULL;

    // every proper prefix must be rejected (and must not read out of bounds)
    int allRejected = 1;
    for (size_t n = 0; n < len; n++) {
        int status = decodeCompressed (buf, n, &d);
        allRejected &= status != SVDC_OK && d == NULL;
        freeCompressed (d);
        d = NULL;
    }
    CHECK (allRejected);

    // trailing garbage
    unsigned char* longer = (unsigned char*) malloc (len + 1);
    memcpy (longer, buf, len);
    longer[len] = 0;
    CHECK (decodeCompressed (longer, len + 1, &d) == SVDC_ERR_CORRUPT);
    free (longer);

    // bad magic, version, bit depth, rank
    unsigned char* copy = (unsigned char*) malloc (len);
    memcpy (copy, buf, len);
    copy[0] = 'X';
    CHECK (decodeCompressed (copy, len, &d) == SVDC_ERR_FORMAT);
    memcpy (copy, buf, len);
    copy[4] = 2;
    CHECK (decodeCompressed (copy, len, &d) == SVDC_ERR_FORMAT);
    memcpy (copy, buf, len);
    copy[5] = 12;
    CHECK (decodeCompressed (copy, len, &d) == SVDC_ERR_CORRUPT);
    memcpy (copy, buf, len);
    copy[17] = 0xFF;  // k of the first channel becomes huge
    CHECK (decodeCompressed (copy, len, &d) != SVDC_OK);

    // flipping any single byte must never crash (result may or may not decode)
    for (size_t i = 0; i < len; i++) {
        memcpy (copy, buf, len);
        copy[i] ^= 0x5A;
        if (decodeCompressed (copy, len, &d) == SVDC_OK) {
            Image* out = decompressImage (d);
            freeImage (out);
        }
        freeCompressed (d);
        d = NULL;
    }
    free (copy);
    free (buf);
    freeCompressed (c);
    freeImage (img);
}

static void testExampleFiles (void) {
    struct {
        const char* name;
        int width, height, channels;
    } files[] = {
        { EXAMPLES_DIR "/test2.png", 735, 526, 3 },  // non-square RGB
        { EXAMPLES_DIR "/test3.png", 360, 360, 4 },  // RGBA
        { EXAMPLES_DIR "/img.jpg", 868, 1541, 3 },   // JPEG, portrait
    };
    for (size_t i = 0; i < sizeof files / sizeof files[0]; i++) {
        Image* img = readImage (files[i].name);
        CHECK (img != NULL);
        if (img == NULL) {
            fprintf (stderr, "  cannot read %s: %s\n", files[i].name, imageFailureReason ());
            continue;
        }
        CHECK (img->width == files[i].width);
        CHECK (img->height == files[i].height);
        CHECK (img->channels == files[i].channels);
        CompressOptions opt = defaultCompressOptions ();
        opt.k = 40;
        Image* out = roundTrip (img, &opt, NULL);
        CHECK (out != NULL && out->width == img->width && out->height == img->height &&
               out->channels == img->channels);
        CHECK (imagePsnr (img, out) > 25.0);
        freeImage (out);
        freeImage (img);
    }
    CHECK (readImage (EXAMPLES_DIR "/does-not-exist.png") == NULL);
}

int main (void) {
    RUN_TEST (testMatrixToChannelClamps);
    RUN_TEST (testFullRankIsLossless);
    RUN_TEST (testQualityGrowsWithK);
    RUN_TEST (testQuantizationLoss);
    RUN_TEST (testConstantAndGrayChannels);
    RUN_TEST (testEnergyMode);
    RUN_TEST (testRejectsBadInput);
    RUN_TEST (testExampleFiles);
    return TEST_RESULT ();
}
