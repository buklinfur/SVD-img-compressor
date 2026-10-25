#include "matrix.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

Matrix* makeMatrix (int rows, int cols) {
    if (rows <= 0 || cols <= 0)
        return NULL;
    if ((size_t) rows > SIZE_MAX / sizeof (double) / (size_t) cols)
        return NULL;

    Matrix* m = (Matrix*) malloc (sizeof (Matrix));
    if (m == NULL)
        return NULL;
    m->data = (double*) calloc ((size_t) rows * (size_t) cols, sizeof (double));
    if (m->data == NULL) {
        free (m);
        return NULL;
    }
    m->rows = rows;
    m->cols = cols;
    return m;
}

void freeMatrix (Matrix* m) {
    if (m == NULL)
        return;
    free (m->data);
    free (m);
}

Matrix* copyMatrix (const Matrix* m) {
    Matrix* res = makeMatrix (m->rows, m->cols);
    if (res == NULL)
        return NULL;
    memcpy (res->data, m->data, matrixSize (m) * sizeof (double));
    return res;
}

Matrix* transposeMatrix (const Matrix* a) {
    Matrix* res = makeMatrix (a->cols, a->rows);
    if (res == NULL)
        return NULL;

    // blocked transpose keeps both reads and writes reasonably cache friendly
    const int block = 32;
    for (int ib = 0; ib < a->rows; ib += block) {
        int iEnd = ib + block < a->rows ? ib + block : a->rows;
        for (int jb = 0; jb < a->cols; jb += block) {
            int jEnd = jb + block < a->cols ? jb + block : a->cols;
            for (int i = ib; i < iEnd; i++)
                for (int j = jb; j < jEnd; j++)
                    matrixSet (res, j, i, matrixGet (a, i, j));
        }
    }
    return res;
}

Matrix* multiply (const Matrix* a, const Matrix* b) {
    if (a->cols != b->rows)
        return NULL;
    Matrix* c = makeMatrix (a->rows, b->cols);
    if (c == NULL)
        return NULL;

    // i-k-j loop order: the inner loop streams over contiguous rows of B and C
    for (int i = 0; i < a->rows; i++) {
        double* cRow = c->data + (size_t) i * (size_t) c->cols;
        const double* aRow = a->data + (size_t) i * (size_t) a->cols;
        for (int k = 0; k < a->cols; k++) {
            double aik = aRow[k];
            if (aik == 0.0)
                continue;
            const double* bRow = b->data + (size_t) k * (size_t) b->cols;
            for (int j = 0; j < b->cols; j++)
                cRow[j] += aik * bRow[j];
        }
    }
    return c;
}

Matrix* multiplyTransA (const Matrix* a, const Matrix* b) {
    if (a->rows != b->rows)
        return NULL;
    Matrix* c = makeMatrix (a->cols, b->cols);
    if (c == NULL)
        return NULL;

    // C[i][j] = sum_k A[k][i] * B[k][j]: accumulate row k of B into row i of C
    for (int k = 0; k < a->rows; k++) {
        const double* aRow = a->data + (size_t) k * (size_t) a->cols;
        const double* bRow = b->data + (size_t) k * (size_t) b->cols;
        for (int i = 0; i < a->cols; i++) {
            double aki = aRow[i];
            if (aki == 0.0)
                continue;
            double* cRow = c->data + (size_t) i * (size_t) c->cols;
            for (int j = 0; j < b->cols; j++)
                cRow[j] += aki * bRow[j];
        }
    }
    return c;
}

Matrix* multiplyTransB (const Matrix* a, const Matrix* b) {
    if (a->cols != b->cols)
        return NULL;
    Matrix* c = makeMatrix (a->rows, b->rows);
    if (c == NULL)
        return NULL;

    // C[i][j] = dot (row i of A, row j of B): both rows are contiguous
    for (int i = 0; i < a->rows; i++) {
        const double* aRow = a->data + (size_t) i * (size_t) a->cols;
        double* cRow = c->data + (size_t) i * (size_t) c->cols;
        for (int j = 0; j < b->rows; j++) {
            const double* bRow = b->data + (size_t) j * (size_t) b->cols;
            double sum = 0.0;
            for (int k = 0; k < a->cols; k++)
                sum += aRow[k] * bRow[k];
            cRow[j] = sum;
        }
    }
    return c;
}

// scaled sum of squares (as in LAPACK dnrm2): returns scale, ssq with sum = scale^2 * ssq
static void sumSquares (const double* x, size_t n, double* scale, double* ssq) {
    for (size_t i = 0; i < n; i++) {
        double v = fabs (x[i]);
        if (v == 0.0)
            continue;
        if (*scale < v) {
            double r = *scale / v;
            *ssq = 1.0 + *ssq * r * r;
            *scale = v;
        } else {
            double r = v / *scale;
            *ssq += r * r;
        }
    }
}

double frobeniusNorm (const Matrix* a) {
    double scale = 0.0, ssq = 1.0;
    sumSquares (a->data, matrixSize (a), &scale, &ssq);
    return scale * sqrt (ssq);
}

double frobeniusNorm2 (const Matrix* a) {
    double sum = 0.0;
    size_t n = matrixSize (a);
    for (size_t i = 0; i < n; i++)
        sum += a->data[i] * a->data[i];
    return sum;
}

double frobeniusDistance (const Matrix* a, const Matrix* b) {
    if (a->rows != b->rows || a->cols != b->cols)
        return (double) NAN;
    double sum = 0.0;
    size_t n = matrixSize (a);
    for (size_t i = 0; i < n; i++) {
        double d = a->data[i] - b->data[i];
        sum += d * d;
    }
    return sqrt (sum);
}

int orthonormalizeColumns (Matrix* a) {
    // work on the transpose so that every column is a contiguous array
    Matrix* t = transposeMatrix (a);
    if (t == NULL)
        return -1;

    const int count = t->rows;
    const size_t len = (size_t) t->cols;
    // a column whose norm drops below this fraction of its original norm after
    // projection is treated as linearly dependent on the previous ones
    const double dependentTol = 1e-10;
    int independent = 0;

    for (int j = 0; j < count; j++) {
        double* x = t->data + (size_t) j * len;
        double original = 0.0;
        for (size_t p = 0; p < len; p++)
            original += x[p] * x[p];
        original = sqrt (original);

        // classical Gram-Schmidt applied twice ("twice is enough", Kahan / Parlett)
        for (int pass = 0; pass < 2; pass++) {
            for (int i = 0; i < j; i++) {
                const double* q = t->data + (size_t) i * len;
                double proj = 0.0;
                for (size_t p = 0; p < len; p++)
                    proj += q[p] * x[p];
                if (proj == 0.0)
                    continue;
                for (size_t p = 0; p < len; p++)
                    x[p] -= proj * q[p];
            }
        }

        double norm = 0.0;
        for (size_t p = 0; p < len; p++)
            norm += x[p] * x[p];
        norm = sqrt (norm);

        if (original == 0.0 || norm <= dependentTol * original) {
            memset (x, 0, len * sizeof (double));
        } else {
            for (size_t p = 0; p < len; p++)
                x[p] /= norm;
            independent++;
        }
    }

    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++)
            matrixSet (a, i, j, matrixGet (t, j, i));
    freeMatrix (t);
    return independent;
}
