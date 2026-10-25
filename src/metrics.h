#ifndef SVDC_METRICS_H
#define SVDC_METRICS_H

#include "image.h"

// mean squared error over all pixels and channels, negative if shapes differ
double imageMse (const Image* a, const Image* b);

// peak signal-to-noise ratio in dB for 8-bit data: 10 log10 (255^2 / MSE).
// INFINITY for identical images, NAN if shapes differ
double imagePsnr (const Image* a, const Image* b);

// wall-clock time in seconds from an arbitrary origin
double wallTime (void);

#endif
