#include "metrics.h"

#include <math.h>
#include <time.h>

double imageMse (const Image* a, const Image* b) {
    if (a->width != b->width || a->height != b->height || a->channels != b->channels)
        return -1.0;
    const size_t count = (size_t) a->width * (size_t) a->height * (size_t) a->channels;
    double sum = 0.0;
    for (size_t i = 0; i < count; i++) {
        double d = (double) a->data[i] - (double) b->data[i];
        sum += d * d;
    }
    return sum / (double) count;
}

double imagePsnr (const Image* a, const Image* b) {
    double mse = imageMse (a, b);
    if (mse < 0.0)
        return (double) NAN;
    if (mse == 0.0)
        return (double) INFINITY;
    return 10.0 * log10 (255.0 * 255.0 / mse);
}

double wallTime (void) {
#if defined(TIME_UTC)
    // C11 timespec_get: wall-clock time with sub-microsecond resolution
    struct timespec ts;
    if (timespec_get (&ts, TIME_UTC) == TIME_UTC)
        return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
#endif
    // fallback for libraries without timespec_get: processor time
    return (double) clock () / (double) CLOCKS_PER_SEC;
}
