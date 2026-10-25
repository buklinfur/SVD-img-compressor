#ifndef SVDC_IMAGE_H
#define SVDC_IMAGE_H

#include "matrix.h"

// 8-bit image with interleaved channels: pixel (x, y), channel c is
// data[(y * width + x) * channels + c]. channels: 1 gray, 2 gray + alpha,
// 3 RGB, 4 RGBA
typedef struct Image {
    int width;
    int height;
    int channels;
    unsigned char* data;
} Image;

// allocate a zero-filled image, NULL on failure or bad size
Image* makeImage (int width, int height, int channels);

// release an image (NULL is allowed)
void freeImage (Image* img);

// load PNG, JPEG, BMP, TGA, GIF (first frame), PSD, PNM. channels are kept as
// stored in the file. returns NULL on failure, see imageFailureReason ()
Image* readImage (const char* filename);

// short description of the last readImage failure
const char* imageFailureReason (void);

// write PNG, JPEG (quality 95), BMP or TGA, chosen by the file extension.
// returns 0 on success
int writeImage (const Image* img, const char* filename);

// 1 if the extension is one writeImage understands
int isWritableImageName (const char* filename);

// copy one channel into a height x width matrix of doubles in [0, 255]
Matrix* channelToMatrix (const Image* img, int channel);

// round, clamp to [0, 255] and store a height x width matrix into one channel
void matrixToChannel (const Matrix* m, Image* img, int channel);

// 1 if every pixel of the channel has the same value (stored in *value)
int isConstantChannel (const Image* img, int channel, unsigned char* value);

#endif
