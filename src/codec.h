#ifndef SVDC_CODEC_H
#define SVDC_CODEC_H

#include <stddef.h>

#include "image.h"
#include "svd.h"

// low-rank representation of one image channel
typedef struct CompressedChannel {
    int constant;         // 1: every pixel equals `value`, no factors stored
    unsigned char value;
    Svd* svd;             // height x width channel ~= U diag (s) V^T (NULL if constant)
} CompressedChannel;

typedef struct Compressed {
    int width;
    int height;
    int channels;
    int bits;  // bits per stored factor value: 8, 16 or 32 (float)
    CompressedChannel channel[4];
} Compressed;

typedef struct CompressOptions {
    int k;           // fixed rank per channel (used when energy <= 0)
    double energy;   // > 0: per channel, the smallest k keeping this fraction of ||A||_F^2
    int bits;        // 8, 16 or 32
    int threads;     // OpenMP threads for the channel loop, 0 = runtime default
    SvdOptions svd;
} CompressOptions;

// k = 50, no energy target, 8-bit factors, randomized SVD
CompressOptions defaultCompressOptions (void);

// error codes shared by the codec functions
enum {
    SVDC_OK = 0,
    SVDC_ERR_ALLOC = -1,
    SVDC_ERR_IO = -2,
    SVDC_ERR_FORMAT = -3,      // not an .svdc file or unsupported version
    SVDC_ERR_TRUNCATED = -4,   // file ends early
    SVDC_ERR_CORRUPT = -5,     // inconsistent header or trailing data
    SVDC_ERR_ARGUMENT = -6
};

const char* codecErrorString (int code);

// factorize every channel. constant channels (e.g. a fully opaque alpha) are
// stored as a single value. returns NULL on allocation failure
Compressed* compressImage (const Image* img, const CompressOptions* opt);

// rebuild the image from the factors, O(width * height * k) per channel
Image* decompressImage (const Compressed* c);

void freeCompressed (Compressed* c);

// serialize to the .svdc byte format (quantizing the factors to c->bits).
// *out is malloc'ed, release it with free ()
int encodeCompressed (const Compressed* c, unsigned char** out, size_t* size);

// parse an .svdc byte buffer. *out receives the dequantized factors
int decodeCompressed (const unsigned char* buf, size_t size, Compressed** out);

int writeCompressedFile (const Compressed* c, const char* filename, size_t* written);
int readCompressedFile (const char* filename, Compressed** out, size_t* size);

// size of the encoded stream in bytes, without building it
size_t encodedSize (const Compressed* c);

#endif
