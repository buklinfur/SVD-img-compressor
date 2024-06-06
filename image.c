#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define STBI_ONLY_PNG
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

static unsigned char zero = '\0';
static int rgbConstructed = 0;
static int rgbInitialized = 0;

typedef struct Image {
    int channels;
    int width;
    int height;
    unsigned char * data;
} Image;

Image* makeImage (int width, int height, int channels, unsigned char* data) {
    Image* img = (Image*) malloc (sizeof(Image) * width * height * channels);
    if (img == NULL) 
        printf ("Failed to allocate Image structure.\n");
    img -> data = realloc (data, sizeof (unsigned char) * width * height * channels);
    if (img -> data == NULL)
        printf ("Failed to allocate Image data.\n");
    img -> height = height;
    img -> width = width;
    img -> channels = channels;
    return img;
}

unsigned char* readImage (const char* filename, int *width, int *height, int *channels) {
    unsigned char* img;
    img = stbi_load(filename, width, height, channels, 3);
    *channels = 3?3:3;
    return img;
}

void printImage (Image* img, const char* fileName) {
    stbi_write_png(fileName, img->width, img->height, img->channels, 
                    img->data, img->width * img->channels);
} 

void makeRgb (Image img, unsigned char* rArr, unsigned char* gArr, unsigned char* bArr) {
    if (rgbConstructed) 
        return;

    int k = 0;
    int rI = 0, gI = 0, bI = 0;
    for (int i = 0; i < img.height; i++) {
        for (int j = 0; j < (img.width * img.channels); j++) {
            if ( k >= 3)
                k = 0;
            if (k == 0) {
                rArr [rI] = img.data[i*img.width*img.channels + j];
                rI++;
            } else if (k == 1) {
                gArr [gI] = img.data[i*img.width*img.channels + j];
                gI++;
            } else if (k == 2) {
                bArr [bI] = img.data[i*img.width*img.channels + j];
                bI++;
            }
            k++;
        }
    }
    rgbConstructed++;
}

void getRgbData (int height, int width, 
                unsigned char* rIn, unsigned char* gIn, unsigned char* bIn,
                unsigned char* rOut, unsigned char* gOut, unsigned char* bOut) {
    if (rgbConstructed == 0) {
        printf ("RGB data not found.\n");
        return;
    }

    int k = 0;
    int rI = 0, gI = 0, bI = 0;

    for (int i = 0; i < height * width * 3; i++) {
            if (k >= 3) 
                k = 0;
            if (k == 0) {
                rOut [i+k] = rIn [rI];
                rI++;
                gOut [i] = zero;
                bOut [i] = zero;
            } else if (k == 1) {
                gOut [i] = gIn [gI];
                gI++;
                rOut [i] = zero;
                bOut [i] = zero;
            } else if (k == 2) {
                bOut [i] = bIn [bI];
                bI++;
                gOut [i] = zero;
                rOut [i] = zero;
            }
            k++;
    }
}

unsigned char* rgbToFull (int width, int height, 
                        unsigned char* r, unsigned char* g, unsigned char* b) {
    if (rgbConstructed == 0) {
        printf ("RGB data not found.\n");
        return NULL;
    }

    int k = 0;
    int rI = 0, gI = 0, bI = 0;

    unsigned char* result = (unsigned char*) malloc (sizeof(unsigned char) * width * height * 3);
    for (int i = 0; i < height * width * 3; i++) {
            if (k >= 3) 
                k = 0;
            if (k == 0) {
                result [i] = r [rI];
                rI++;
            } else if (k == 1) {
                result [i] = g [gI];
                gI++;
            } else if (k == 2) {
                result [i] = b [bI];
                bI++;
            }
            k++;
    }

    return result;
}