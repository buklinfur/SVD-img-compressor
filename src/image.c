#include "image.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define STBI_NO_HDR
#define STBI_NO_LINEAR
#include "stb_image.h"
#include "stb_image_write.h"

// stb_image refuses anything larger, keep the same limit everywhere
#define SVDC_MAX_SIDE (1 << 24)

Image* makeImage (int width, int height, int channels) {
    if (width <= 0 || height <= 0 || width > SVDC_MAX_SIDE || height > SVDC_MAX_SIDE)
        return NULL;
    if (channels < 1 || channels > 4)
        return NULL;
    if ((size_t) width > SIZE_MAX / (size_t) height / (size_t) channels)
        return NULL;

    Image* img = (Image*) malloc (sizeof (Image));
    if (img == NULL)
        return NULL;
    img->data = (unsigned char*) calloc ((size_t) width * (size_t) height * (size_t) channels, 1);
    if (img->data == NULL) {
        free (img);
        return NULL;
    }
    img->width = width;
    img->height = height;
    img->channels = channels;
    return img;
}

void freeImage (Image* img) {
    if (img == NULL)
        return;
    free (img->data);
    free (img);
}

Image* readImage (const char* filename) {
    int width = 0, height = 0, channels = 0;
    // desired channels = 0: keep what the file stores (1..4)
    unsigned char* pixels = stbi_load (filename, &width, &height, &channels, 0);
    if (pixels == NULL)
        return NULL;

    // copy into our own allocation so that every Image is released with freeImage
    Image* img = makeImage (width, height, channels);
    if (img != NULL)
        memcpy (img->data, pixels, (size_t) width * (size_t) height * (size_t) channels);
    stbi_image_free (pixels);
    return img;
}

const char* imageFailureReason (void) {
    const char* reason = stbi_failure_reason ();
    return reason != NULL ? reason : "unknown error";
}

// lower-case extension of a file name ("png", "jpg"...), empty if none
static void extensionOf (const char* filename, char* ext, size_t size) {
    ext[0] = '\0';
    const char* dot = strrchr (filename, '.');
    const char* slash = strrchr (filename, '/');
    const char* backslash = strrchr (filename, '\\');
    if (dot == NULL || (slash != NULL && dot < slash) || (backslash != NULL && dot < backslash))
        return;
    size_t i = 0;
    for (const char* p = dot + 1; *p != '\0' && i + 1 < size; p++, i++)
        ext[i] = (char) tolower ((unsigned char) *p);
    ext[i] = '\0';
}

int isWritableImageName (const char* filename) {
    char ext[8];
    extensionOf (filename, ext, sizeof ext);
    return strcmp (ext, "png") == 0 || strcmp (ext, "jpg") == 0 || strcmp (ext, "jpeg") == 0 ||
           strcmp (ext, "bmp") == 0 || strcmp (ext, "tga") == 0;
}

int writeImage (const Image* img, const char* filename) {
    char ext[8];
    extensionOf (filename, ext, sizeof ext);
    int ok = 0;
    if (strcmp (ext, "png") == 0)
        ok = stbi_write_png (filename, img->width, img->height, img->channels, img->data,
                             img->width * img->channels);
    else if (strcmp (ext, "jpg") == 0 || strcmp (ext, "jpeg") == 0)
        ok = stbi_write_jpg (filename, img->width, img->height, img->channels, img->data, 95);
    else if (strcmp (ext, "bmp") == 0)
        ok = stbi_write_bmp (filename, img->width, img->height, img->channels, img->data);
    else if (strcmp (ext, "tga") == 0)
        ok = stbi_write_tga (filename, img->width, img->height, img->channels, img->data);
    return ok ? 0 : -1;
}

Matrix* channelToMatrix (const Image* img, int channel) {
    Matrix* m = makeMatrix (img->height, img->width);
    if (m == NULL)
        return NULL;
    const size_t stride = (size_t) img->channels;
    const size_t count = matrixSize (m);
    for (size_t i = 0; i < count; i++)
        m->data[i] = (double) img->data[i * stride + (size_t) channel];
    return m;
}

void matrixToChannel (const Matrix* m, Image* img, int channel) {
    const size_t stride = (size_t) img->channels;
    const size_t count = matrixSize (m);
    for (size_t i = 0; i < count; i++) {
        double v = floor (m->data[i] + 0.5);
        if (!(v > 0.0))  // also maps NaN to 0
            v = 0.0;
        else if (v > 255.0)
            v = 255.0;
        img->data[i * stride + (size_t) channel] = (unsigned char) v;
    }
}

int isConstantChannel (const Image* img, int channel, unsigned char* value) {
    const size_t stride = (size_t) img->channels;
    const size_t count = (size_t) img->width * (size_t) img->height;
    const unsigned char first = img->data[(size_t) channel];
    for (size_t i = 1; i < count; i++)
        if (img->data[i * stride + (size_t) channel] != first)
            return 0;
    if (value != NULL)
        *value = first;
    return 1;
}
