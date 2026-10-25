// accuracy tests for the Jacobi and randomized SVD on matrices with known answers

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "matrix.h"
#include "svd.h"
#include "test_util.h"

static uint64_t rngState = 12345;

static double randomUniform (void) {
    uint64_t z = (rngState += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return (double) (z >> 11) / 9007199254740992.0 * 2.0 - 1.0;  // [-1, 1)
}

static Matrix* randomMatrix (int rows, int cols) {
    Matrix* m = makeMatrix (rows, cols);
    for (size_t i = 0; i < matrixSize (m); i++)
        m->data[i] = randomUniform ();
    return m;
}

// rows x cols matrix with orthonormal columns (cols <= rows)
static Matrix* randomOrthonormal (int rows, int cols) {
    Matrix* m = randomMatrix (rows, cols);
    orthonormalizeColumns (m);
    return m;
}

// A = Q1 * diag (s) * Q2^T with random orthonormal Q1 (m x r), Q2 (n x r)
static Matrix* matrixWithSpectrum (int m, int n, const double* s, int r) {
    Matrix* q1 = randomOrthonormal (m, r);
    Matrix* q2 = randomOrthonormal (n, r);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < r; j++)
            matrixSet (q1, i, j, matrixGet (q1, i, j) * s[j]);
    Matrix* a = multiplyTransB (q1, q2);
    freeMatrix (q1);
    freeMatrix (q2);
    return a;
}

// max |Q^T Q - I| over columns that belong to nonzero singular values
static double orthogonalityError (const Matrix* q, const double* s, int rank) {
    double worst = 0.0;
    for (int a = 0; a < rank; a++) {
        if (s[a] == 0.0)
            continue;
        for (int b = 0; b < rank; b++) {
            if (s[b] == 0.0)
                continue;
            double d = 0.0;
            for (int i = 0; i < q->rows; i++)
                d += matrixGet (q, i, a) * matrixGet (q, i, b);
            double err = fabs (d - (a == b ? 1.0 : 0.0));
            if (err > worst)
                worst = err;
        }
    }
    return worst;
}

static void checkDecomposition (const Matrix* a, const Svd* svd, double tol) {
    CHECK (svd != NULL);
    if (svd == NULL)
        return;
    for (int i = 0; i < svd->rank; i++) {
        CHECK (svd->s[i] >= 0.0);
        if (i > 0)
            CHECK (svd->s[i] <= svd->s[i - 1]);
    }
    Matrix* r = svdReconstruct (svd, 0);
    double norm = frobeniusNorm (a);
    double err = frobeniusDistance (a, r);
    CHECK_LESS (norm > 0.0 ? err / norm : err, tol);
    CHECK_LESS (orthogonalityError (svd->u, svd->s, svd->rank), tol);
    CHECK_LESS (orthogonalityError (svd->v, svd->s, svd->rank), tol);
    freeMatrix (r);
}

static void testDiagonal (void) {
    const double diag[5] = { 3.0, -7.0, 0.5, 10.0, 2.0 };
    const double expected[5] = { 10.0, 7.0, 3.0, 2.0, 0.5 };
    Matrix* a = makeMatrix (5, 5);
    for (int i = 0; i < 5; i++)
        matrixSet (a, i, i, diag[i]);
    Svd* svd = svdJacobi (a, 0, 0, NULL);
    checkDecomposition (a, svd, 1e-14);
    for (int i = 0; i < 5; i++)
        CHECK_NEAR (svd->s[i], expected[i], 1e-14);
    freeSvd (svd);
    freeMatrix (a);
}

static void testRectangularDiagonal (void) {
    // 6 x 4 and its transpose: singular values are |diagonal|, sorted
    Matrix* a = makeMatrix (6, 4);
    matrixSet (a, 0, 0, 1.0);
    matrixSet (a, 1, 1, -4.0);
    matrixSet (a, 2, 2, 2.0);
    matrixSet (a, 3, 3, 3.0);
    Matrix* at = transposeMatrix (a);
    const double expected[4] = { 4.0, 3.0, 2.0, 1.0 };

    Svd* s1 = svdJacobi (a, 0, 0, NULL);
    Svd* s2 = svdJacobi (at, 0, 0, NULL);
    checkDecomposition (a, s1, 1e-14);
    checkDecomposition (at, s2, 1e-14);
    CHECK (s1->rank == 4 && s2->rank == 4);
    CHECK (s1->u->rows == 6 && s1->v->rows == 4);
    CHECK (s2->u->rows == 4 && s2->v->rows == 6);
    for (int i = 0; i < 4; i++) {
        CHECK_NEAR (s1->s[i], expected[i], 1e-14);
        CHECK_NEAR (s2->s[i], expected[i], 1e-14);
    }
    freeSvd (s1);
    freeSvd (s2);
    freeMatrix (a);
    freeMatrix (at);
}

static void testOrthogonal (void) {
    // every singular value of an orthogonal matrix equals 1
    Matrix* q = randomOrthonormal (12, 12);
    Svd* svd = svdJacobi (q, 0, 0, NULL);
    checkDecomposition (q, svd, 1e-13);
    for (int i = 0; i < svd->rank; i++)
        CHECK_NEAR (svd->s[i], 1.0, 1e-13);
    freeSvd (svd);
    freeMatrix (q);
}

static void testRankOne (void) {
    // A = x y^T has the single nonzero singular value ||x|| ||y||
    const int m = 9, n = 7;
    Matrix* x = randomMatrix (m, 1);
    Matrix* y = randomMatrix (n, 1);
    Matrix* a = multiplyTransB (x, y);
    double expected = frobeniusNorm (x) * frobeniusNorm (y);

    Svd* svd = svdJacobi (a, 0, 0, NULL);
    checkDecomposition (a, svd, 1e-13);
    CHECK_NEAR (svd->s[0], expected, 1e-13 * expected);
    for (int i = 1; i < svd->rank; i++)
        CHECK_LESS (svd->s[i], 1e-13 * expected);

    Svd* rnd = svdRandomized (a, 1, NULL);
    CHECK (rnd != NULL && rnd->rank == 1);
    CHECK_NEAR (rnd->s[0], expected, 1e-13 * expected);

    freeSvd (svd);
    freeSvd (rnd);
    freeMatrix (a);
    freeMatrix (x);
    freeMatrix (y);
}

static void testZero (void) {
    Matrix* a = makeMatrix (3, 5);
    Svd* svd = svdJacobi (a, 0, 0, NULL);
    CHECK (svd != NULL && svd->rank == 3);
    for (int i = 0; i < svd->rank; i++)
        CHECK (svd->s[i] == 0.0);
    Matrix* r = svdReconstruct (svd, 0);
    CHECK (frobeniusNorm (r) == 0.0);
    freeMatrix (r);
    freeSvd (svd);

    Svd* e = svdForEnergy (a, 0.95, NULL);
    CHECK (e != NULL && e->rank == 1 && e->s[0] == 0.0);
    freeSvd (e);
    freeMatrix (a);
}

static void testRandomShapes (void) {
    const int shapes[][2] = { { 40, 25 }, { 25, 40 }, { 33, 33 }, { 1, 17 }, { 17, 1 }, { 2, 2 }, { 120, 7 } };
    for (size_t t = 0; t < sizeof (shapes) / sizeof (shapes[0]); t++) {
        Matrix* a = randomMatrix (shapes[t][0], shapes[t][1]);
        JacobiStats stats;
        Svd* svd = svdJacobi (a, 0, 0, &stats);
        CHECK (stats.converged);
        checkDecomposition (a, svd, 1e-12);
        freeSvd (svd);
        freeMatrix (a);
    }
}

static void testKnownSpectrum (void) {
    // widely spread spectrum: sigma_i = 10^(-i/2), i = 0..23 (condition number ~1e11)
    double s[24];
    for (int i = 0; i < 24; i++)
        s[i] = pow (10.0, -0.5 * i);
    Matrix* a = matrixWithSpectrum (60, 45, s, 24);
    Svd* svd = svdJacobi (a, 0, 0, NULL);
    checkDecomposition (a, svd, 1e-12);
    for (int i = 0; i < 24; i++)
        CHECK_NEAR (svd->s[i], s[i], 1e-14);
    for (int i = 24; i < svd->rank; i++)
        CHECK_LESS (svd->s[i], 1e-14);
    freeSvd (svd);
    freeMatrix (a);
}

static void testTruncation (void) {
    Matrix* a = randomMatrix (30, 20);
    Svd* full = svdJacobi (a, 0, 0, NULL);
    Svd* top = svdJacobi (a, 5, 0, NULL);
    CHECK (top->rank == 5 && top->u->cols == 5 && top->v->cols == 5);
    for (int i = 0; i < 5; i++)
        CHECK_NEAR (top->s[i], full->s[i], 1e-14);

    // Eckart-Young: ||A - A_k||_F^2 = sum of the discarded sigma_i^2
    double tail = 0.0;
    for (int i = 5; i < full->rank; i++)
        tail += full->s[i] * full->s[i];
    Matrix* r = svdReconstruct (top, 0);
    CHECK_NEAR (frobeniusDistance (a, r), sqrt (tail), 1e-12);

    truncateSvd (full, 5);
    CHECK (full->rank == 5);
    Matrix* r2 = svdReconstruct (full, 0);
    CHECK_LESS (frobeniusDistance (r, r2), 1e-12);

    freeMatrix (r);
    freeMatrix (r2);
    freeSvd (full);
    freeSvd (top);
    freeMatrix (a);
}

static void testRandomizedDecaying (void) {
    // geometric decay similar to natural images
    const int m = 300, n = 200, r = 200;
    double s[200];
    for (int i = 0; i < r; i++)
        s[i] = 100.0 * pow (0.85, i);
    Matrix* a = matrixWithSpectrum (m, n, s, r);

    const int k = 20;
    SvdOptions opt = defaultSvdOptions ();
    Svd* svd = svdRandomized (a, k, &opt);
    CHECK (svd != NULL && svd->rank == k);
    checkDecomposition (a, svd, 1.0);  // ordering and sign checks only

    for (int i = 0; i < k; i++)
        CHECK_NEAR (svd->s[i], s[i], 1e-8 * s[0]);

    double optimal = 0.0;
    for (int i = k; i < r; i++)
        optimal += s[i] * s[i];
    optimal = sqrt (optimal);
    Matrix* rec = svdReconstruct (svd, 0);
    double err = frobeniusDistance (a, rec);
    CHECK_LESS (err, 1.0001 * optimal);
    CHECK_LESS (orthogonalityError (svd->u, svd->s, svd->rank), 1e-12);

    freeMatrix (rec);
    freeSvd (svd);
    freeMatrix (a);
}

static void testRandomizedFallback (void) {
    // k + oversampling >= min (m, n): must behave like the exact SVD
    Matrix* a = randomMatrix (15, 12);
    Svd* exact = svdJacobi (a, 0, 0, NULL);
    Svd* rnd = svdRandomized (a, 8, NULL);
    CHECK (rnd->rank == 8);
    for (int i = 0; i < 8; i++)
        CHECK_NEAR (rnd->s[i], exact->s[i], 1e-13);
    freeSvd (exact);
    freeSvd (rnd);
    freeMatrix (a);
}

static void testEnergyRank (void) {
    const double s[4] = { 3.0, 2.0, 1.0, 1.0 };  // energies 9, 4, 1, 1 of 15
    CHECK (rankForEnergy (s, 4, 15.0, 0.5) == 1);
    CHECK (rankForEnergy (s, 4, 15.0, 0.61) == 2);
    CHECK (rankForEnergy (s, 4, 15.0, 13.0 / 15.0) == 2);
    CHECK (rankForEnergy (s, 4, 15.0, 0.9) == 3);
    CHECK (rankForEnergy (s, 4, 15.0, 1.0) == 4);

    double spec[60];
    for (int i = 0; i < 60; i++)
        spec[i] = pow (0.9, i);
    Matrix* a = matrixWithSpectrum (200, 150, spec, 60);
    SvdOptions jac = defaultSvdOptions ();
    jac.method = SVD_JACOBI;
    Svd* e1 = svdForEnergy (a, 0.99, &jac);
    Svd* e2 = svdForEnergy (a, 0.99, NULL);
    int expected = rankForEnergy (spec, 60, frobeniusNorm2 (a), 0.99);
    CHECK (e1->rank == expected);
    CHECK (e2->rank == expected);
    freeSvd (e1);
    freeSvd (e2);
    freeMatrix (a);
}

int main (void) {
    RUN_TEST (testDiagonal);
    RUN_TEST (testRectangularDiagonal);
    RUN_TEST (testOrthogonal);
    RUN_TEST (testRankOne);
    RUN_TEST (testZero);
    RUN_TEST (testRandomShapes);
    RUN_TEST (testKnownSpectrum);
    RUN_TEST (testTruncation);
    RUN_TEST (testRandomizedDecaying);
    RUN_TEST (testRandomizedFallback);
    RUN_TEST (testEnergyRank);
    return TEST_RESULT ();
}
