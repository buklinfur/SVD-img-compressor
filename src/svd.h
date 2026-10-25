#ifndef SVDC_SVD_H
#define SVDC_SVD_H

#include <stdint.h>

#include "matrix.h"

// thin (truncated) singular value decomposition A ~= U * diag (s) * V^T
//   u: rows x rank, orthonormal columns (columns for zero singular values are zero)
//   s: rank singular values, sorted in non-increasing order, all >= 0
//   v: cols x rank, orthonormal columns (V itself, not V^T)
typedef struct Svd {
    int rows;
    int cols;
    int rank;
    Matrix* u;
    double* s;
    Matrix* v;
} Svd;

typedef enum SvdMethod {
    SVD_JACOBI = 0,     // exact one-sided Jacobi (Hestenes), O(m n^2) per sweep
    SVD_RANDOMIZED = 1  // randomized range finder + Jacobi, O(m n k)
} SvdMethod;

typedef struct JacobiStats {
    int sweeps;      // sweeps performed
    long rotations;  // plane rotations applied in total
    int converged;   // 1 if the stopping criterion was met before maxSweeps
} JacobiStats;

typedef struct SvdOptions {
    SvdMethod method;
    int oversampling;  // randomized: extra columns in the sketch (default 10)
    int powerIters;    // randomized: power iterations (default 2)
    uint64_t seed;     // randomized: seed of the Gaussian test matrix
    int maxSweeps;     // Jacobi: upper bound on sweeps (default 40)
} SvdOptions;

// default options: randomized, p = 10, q = 2, seed = 1, maxSweeps = 40
SvdOptions defaultSvdOptions (void);

// exact thin SVD with one-sided Jacobi rotations.
// k <= 0 or k > min (rows, cols) keeps all min (rows, cols) triplets.
// stats may be NULL. returns NULL on allocation failure.
Svd* svdJacobi (const Matrix* a, int k, int maxSweeps, JacobiStats* stats);

// approximate top-k SVD (Halko, Martinsson, Tropp 2011, Alg. 4.4 + 5.1).
// falls back to the exact Jacobi SVD when k + oversampling covers the full rank.
Svd* svdRandomized (const Matrix* a, int k, const SvdOptions* opt);

// dispatch on opt->method (opt may be NULL for defaults)
Svd* svdCompute (const Matrix* a, int k, const SvdOptions* opt);

// smallest k such that s[0..k-1] holds at least `fraction` of `totalEnergy`
// (energy = sum of squared singular values = ||A||_F^2). returns count if never reached.
int rankForEnergy (const double* s, int count, double totalEnergy, double fraction);

// top-k SVD whose k is the smallest rank capturing `fraction` of ||A||_F^2.
// for the randomized method the sketch size grows until the target is met.
Svd* svdForEnergy (const Matrix* a, double fraction, const SvdOptions* opt);

// keep only the first k triplets (in place). k >= rank is a no-op
int truncateSvd (Svd* svd, int k);

// U_k * diag (s_k) * V_k^T for the first k triplets (k <= 0 means all), O(m n k)
Matrix* svdReconstruct (const Svd* svd, int k);

void freeSvd (Svd* svd);

#endif
