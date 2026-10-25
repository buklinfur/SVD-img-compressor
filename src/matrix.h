#ifndef SVDC_MATRIX_H
#define SVDC_MATRIX_H

#include <stddef.h>

// dense row-major matrix of doubles: element (i, j) is data[i * cols + j]
typedef struct Matrix {
    int rows;
    int cols;
    double* data;
} Matrix;

// allocate a zero-filled rows x cols matrix, NULL on failure or bad size
Matrix* makeMatrix (int rows, int cols);

// release a matrix (NULL is allowed)
void freeMatrix (Matrix* m);

// deep copy, NULL on failure
Matrix* copyMatrix (const Matrix* m);

// element access helpers
static inline double matrixGet (const Matrix* m, int i, int j) {
    return m->data[(size_t) i * (size_t) m->cols + (size_t) j];
}

static inline void matrixSet (Matrix* m, int i, int j, double value) {
    m->data[(size_t) i * (size_t) m->cols + (size_t) j] = value;
}

// number of elements
static inline size_t matrixSize (const Matrix* m) {
    return (size_t) m->rows * (size_t) m->cols;
}

// B = A^T
Matrix* transposeMatrix (const Matrix* a);

// C = A * B
Matrix* multiply (const Matrix* a, const Matrix* b);

// C = A^T * B (A is k x m, B is k x n, C is m x n)
Matrix* multiplyTransA (const Matrix* a, const Matrix* b);

// C = A * B^T (A is m x k, B is n x k, C is m x n)
Matrix* multiplyTransB (const Matrix* a, const Matrix* b);

// Frobenius norm, computed with scaling to avoid overflow
double frobeniusNorm (const Matrix* a);

// squared Frobenius norm (plain sum of squares)
double frobeniusNorm2 (const Matrix* a);

// ||A - B||_F, matrices must have equal shape
double frobeniusDistance (const Matrix* a, const Matrix* b);

// orthonormalize the columns of A in place (Gram-Schmidt applied twice).
// columns that turn out numerically dependent on the previous ones are set to zero.
// returns the number of nonzero (independent) columns, or -1 on allocation failure
int orthonormalizeColumns (Matrix* a);

#endif
