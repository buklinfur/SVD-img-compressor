#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct Image {
    int channels;
    int width;
    int height;
    unsigned char * data;
} Image;

// data allocation and assembly into the Image structure
Image* makeImage (int width, int height, int channels, unsigned char* data);

// image input
unsigned char* readImage (const char* filename, int *width, int *height, int *channels);

// image output
void printImage (Image* img, const char* fileName);

// creating distinct rgb mono arrays
void makeRgb (Image img, unsigned char* rArr, unsigned char* gArr, unsigned char* bArr);

// converting rgb mono arrays into proper mono image data
void getRgbData (int height, int width, 
                unsigned char* rIn, unsigned char* gIn, unsigned char* bIn,
                unsigned char* rOut, unsigned char* gOut, unsigned char* bOut);

// converting rgb mono arrays into proper image data
unsigned char* rgbToFull (int width, int height, 
                        unsigned char* r, unsigned char* g, unsigned char* b);