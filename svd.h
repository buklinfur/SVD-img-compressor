#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef struct Matrix {
    unsigned char* data;
    int width; // m
    int height; // n
} Matrix;

// Matrix memory allocation
Matrix* makeMatrix (int m, int n, int dataSize);

// matrix construction from data array
void dataToMatrix (Matrix* matrix, unsigned char* data, int n, int m, int dataSize);

// matrix output (mostly for debugging)
void print_matrix(FILE * f, Matrix* m);

// get diagonal matrix
double* diag (double* A, int row_b, int column_b);

// get transposed matrix
void tran(double* A, int row, int column);

// get U, S, V matrices ov SVD
void svd(double* A, double* U, double* S, double* V, int row, int column);

// restore image data from U*W*VT matrices, where k - a number of singular values
// the higher the number - the better the quality
Matrix* svdCompression (Matrix* A, int k); 