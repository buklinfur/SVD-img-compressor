#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <memory.h>

#include <lapacke.h>
#include <lapack.h>
#include <cblas.h>

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

typedef struct Matrix
{
    unsigned char *data;
    int width;
    int height;
} Matrix;

Matrix *makeMatrix(int width, int height, int dataSize)
{
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix) * width * height * dataSize);
    if (matrix == NULL)
        printf("Failed to allocate Matrix structure");
    matrix->data = (unsigned char *)malloc(sizeof(unsigned char) * height * width * dataSize);
    matrix->height = height;
    matrix->width = width * dataSize;
    return matrix;
}

void dataToMatrix(Matrix *matrix, unsigned char *data, int height,
                  int width, int dataSize)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < (width * dataSize); j++)
        {
            matrix->data[i * width * dataSize + j] = data[i * width * dataSize + j];
        }
    }
}

void print_matrix(FILE *f, Matrix *m)
{
    for (int i = 0; i < m->height; i++)
    {
        for (int j = 0; j < (m->width); j++)
        {
            fwrite(&m->data[i * m->width + j], sizeof(unsigned char), 1, f);
        }
    }
}

double* toDouble (unsigned char* data, int n, int m) {
    double* res = (double*)malloc(sizeof(double) * n * m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res [m*i + j] = (double) data [m*i + j];
        }
    } 
    return res;
}

unsigned char* toChar (double* data, int n, int m) {
    unsigned char* res = (unsigned char*) malloc (sizeof(unsigned char) * n * m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double value = round(data [m*i + j]);
            if (value < 0) {
                if (abs(value)*2 > 255) 
                    value = abs(value);
                else
                    value = 0;
            } else if (value > 255) {
                value = 255;
            }
            res [m*i + j] = (unsigned char) value;
        }
    } 
    return res;
}

double* diag (double* A, int row_b, int column_b) {
	double* B = (double*) malloc (sizeof(double)*row_b*column_b*column_b);
	for (int i = 0; i < row_b; i++) {
		for (int j = 0; j < column_b; j++) {
			if (j == i) 
				B[i*column_b+j] = A[i];
		}
	}
    return B;
}

double* tran (double* A, int row, int column) {
	double* B = (double*) malloc (sizeof(double)*row*column);
	double* transpose;
	double* ptr_A = A;

	for (int i = 0; i < row; i++) {
		transpose = &B[i];
		for (int j = 0; j < column; j++) {
			*transpose = *ptr_A;
			ptr_A++;
			transpose += row;
		}
	}
	return B;
}

/*
sizes:
U = row*k
S = k*k
V = k*column

need to implement RES = U*S*V, then:
res ij = sum [r=1..k] (sum [y=1..k] (uik*syr) * vrj), where:
*/

/*
double* getMult (double* U, double* S, double* V, int row, int column, int k) {
    double* res = (double*) malloc (sizeof(double)*row*column);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            res [i * column + j]  = 0;
            for (int r = 0; r < k; r++) {
                int intermediate = 0;
                for (int y = 0; y < k; y++) {
                    intermediate += U [i * column + y] * S [y * row + r];
                }
                res [i * column + j] += intermediate * V [r * row + k];
            }
        }
    }
    return res;
}
*/

double* getMult (double* U, double* S, double* V, int row, int column) {
    double* intermediate = (double*) malloc (sizeof(double)*row*column);
    double* res = (double*) malloc (sizeof(double)*row*column);

    const int R = row;
    const int K = fmin (row, column);
    const int C = column;
    const double alpha = 1;
    const double beta = 0;

    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, R, K, K, alpha, U, R, S, K, beta, intermediate, R);
    cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, R, C, K, alpha, intermediate, R, V, K, beta, res, R);

    return res;
}

double* shrinker (double* A, int rows, int cols, int rowsNew, int colsNew) {
    double* res = (double*) malloc (sizeof(double) * rowsNew * colsNew);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (i<rowsNew && j<colsNew) 
                res [colsNew*i + j] = A [cols*i + j];
            else 
                res [cols*i + j] = 0;
        }
    }
    return res;
}

double* svd(double* A, double* U, double* S, double* V, int row, int column) {
	lapack_int M = row;
	lapack_int N = column;
    double* superb;

 	printf ("\t\tAttempting to use LAPACKE function...\n");
	lapack_int i = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A', 'A', M, N, A, N, S, U, N, V, N, superb);
    printf ("\t\tLAPACKE returned: %d\n", i);
 	printf ("\t\tLAPACKE function complete.\n");

 	printf ("\t\tAttempting to make S diagonal...\n");
	double* S_ = diag(S, row, column);
 	printf ("\t\tDiagonalization complete.\n");
    free (S);
    return S_;
}

Matrix* svdCompression (Matrix* A, int k) {
    int column = A->width;
    int row = A->height;

    printf ("\tAttempting to allocate matrices for SVD...\n");
    double* matrix = toDouble (A->data, row, column);
    double* U = (double*) malloc (sizeof(double)*row*row);
    double* S = (double*) malloc (sizeof(double)*fmin(row, column));
    double* V = (double*) malloc (sizeof(double)*column*column);
    printf ("\tAllocation complete.\n");

    printf ("\tAttempting to apply SVD to the input matrix...\n");
    double* S_ = shrinker(svd (matrix, U, S, V, row, column), row, column, k, k);
    printf ("\tApplying SVD to the input matrix complete.\n");
    
    printf ("\tAttempting to multiply singular matrices...\n");
    double* resArrDouble = getMult (U, S_, V, row, column);
    unsigned char* resArrChar = toChar (resArrDouble, row, column);
    printf ("\tMultiplication complete.\n");

    Matrix* res = makeMatrix (column, row, 1);
    dataToMatrix (res, resArrChar, row, column, 1);

    return res;
}

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

int main () {
    FILE* hexout = fopen ("output.out", "wb"); // for debuggin image data
    FILE* output = fopen ("output.txt", "w");

    //inputValidation (argc);

    int standardCompression = 100;
    //int* compressionValue = (argc == 2) ? &standardCompression : ((int*) argv[2]);
    int compressionValue = standardCompression;

// image load sequence    
    int width, height, channels;
    unsigned char *data = readImage ("test.png", &width, &height, &channels);
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

    printf ("Attempting split rgb image creation...\n");
    unsigned char *rIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    unsigned char *gIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    unsigned char *bIm = (unsigned char*) malloc (sizeof(unsigned char) * height * width * channels);
    getRgbData (height, width, r, g, b, rIm, gIm, bIm);

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

    free (r);
    free (g);
    free (b);

    printf ("Attempting to apply SVD compression...\n");
    Matrix* resR = svdCompression (matrixR, compressionValue);
    Matrix* resG = svdCompression (matrixG, compressionValue);
    Matrix* resB = svdCompression (matrixB, compressionValue);
    printf ("SVD compression complete\n\n");

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