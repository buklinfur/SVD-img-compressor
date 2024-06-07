#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "image.h"
#include "svd.h"

void inputValidation (int argc) {
    if(argc == 1) {
        printf("No image input.\n\
Insert the filename as the command line argument.\n");
        exit(6);
    }
    if(argc == 2) {
        printf("No compression parameter input.\n\
Standard parameter used: 20.\n\
If you want to change it, please consider putting it as the command line argument after filename.\n\n");
    }
}

int main (int argc, char **argv) {
    //FILE* hexout = fopen ("output.out", "wb"); // for debuggin image data
    //FILE* output = fopen ("output.txt", "w");

    inputValidation (argc);

    int standardCompression = 20;
    int compressionValue = argc == 2 ? standardCompression :  atoi (argv[2]);
    

    int* parameter = (int*) argv[2];
// image load sequence    
    int width, height, channels;
    unsigned char *data = readImage (argv[1], &width, &height, &channels);
// image load validation
    if(data == NULL) {
        printf("Error in loading the image\n");
        exit(1);
    }
// data assembly into a structure for a better accessibility
    Image* img = makeImage (width, height, channels, data); 
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n\n", img->width, img->height, img->channels);

    printf ("Processing rgb split...\n");
// rgb mono arrays initialization
    unsigned char *r = (unsigned char*) malloc (sizeof(unsigned char) * height * width);
    unsigned char *g = (unsigned char*) malloc (sizeof(unsigned char) * height * width);
    unsigned char *b = (unsigned char*) malloc (sizeof(unsigned char) * height * width);
    makeRgb (*img, r, g, b);
    printf ("Rgb split complete.\n\n");

    printf ("Attempting to put it back together...\n");
    Image* fullRes = makeImage (width, height, channels, rgbToFull (width, height, r, g, b));
    printf ("Merge complete.\n");
    printImage (fullRes, "testFull.png");
    printf ("Merged image output complete.\n\n");

    printf ("Making matrices for split rgb data...\n");
    Matrix* matrixR = makeMatrix (width, height, 1);
    dataToMatrix (matrixR, r, height, width, 1);
    Matrix* matrixG = makeMatrix (width, height, 1);
    dataToMatrix (matrixG, g, height, width, 1);
    Matrix* matrixB = makeMatrix (width, height, 1);
    dataToMatrix (matrixB, b, height, width, 1);
    printf ("Making matrices complete.\n\n");

    printf ("Attempting to apply SVD compression...\n");
    Matrix* resR = svdCompression (matrixR, compressionValue);
    Matrix* resG = svdCompression (matrixG, compressionValue);
    Matrix* resB = svdCompression (matrixB, compressionValue);
    printf ("SVD compression complete\n\n");

    printf ("Attempting split rgb image creation...\n");
    unsigned char *rIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    unsigned char *gIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    unsigned char *bIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    getRgbData (height, width, resR->data, resG->data, resB->data, rIm, gIm, bIm);

    // rgb nono images initialization
    Image* redImg = makeImage (width, height, 3, rIm);
    Image* greenImg = makeImage (width, height, 3, gIm);
    Image* blueImg = makeImage (width, height, 3, bIm);
    printf ("Split rgb image data conversion complete.\n");

    // rgb mono matrices out (for debugging)
    printImage (redImg, "testRed.png");
    printImage (greenImg, "testGreen.png");
    printImage (blueImg, "testBlue.png");
    printf ("Split rgb image output complete.\n\n");

    printf ("Merging result image...\n");
    unsigned char* res = rgbToFull (height, width, resR->data, resG->data, resB->data);
    printf ("Merging complete.\n\n");

    free (resR->data);
    free (resG->data);
    free (resB->data);
    
    Image* outImage = makeImage (resR->width, resR->height, 3, res);
    printImage (outImage, "testOut.png");
    printf ("Result output complete.\n");

    return 0;
}