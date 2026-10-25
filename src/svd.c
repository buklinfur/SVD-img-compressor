#include "svd.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SVDC_TWO_PI 6.283185307179586476925286766559

SvdOptions defaultSvdOptions (void) {
    SvdOptions opt;
    opt.method = SVD_RANDOMIZED;
    opt.oversampling = 10;
    opt.powerIters = 2;
    opt.seed = 1;
    opt.maxSweeps = 40;
    return opt;
}

void freeSvd (Svd* svd) {
    if (svd == NULL)
        return;
    freeMatrix (svd->u);
    freeMatrix (svd->v);
    free (svd->s);
    free (svd);
}

Svd* makeSvd (int rows, int cols, int rank) {
    if (rank <= 0)
        return NULL;
    Svd* svd = (Svd*) calloc (1, sizeof (Svd));
    if (svd == NULL)
        return NULL;
    svd->rows = rows;
    svd->cols = cols;
    svd->rank = rank;
    svd->u = makeMatrix (rows, rank);
    svd->v = makeMatrix (cols, rank);
    svd->s = (double*) calloc ((size_t) rank, sizeof (double));
    if (svd->u == NULL || svd->v == NULL || svd->s == NULL) {
        freeSvd (svd);
        return NULL;
    }
    return svd;
}

// four independent partial sums: lets the compiler vectorize the reduction
// without -ffast-math and is slightly more accurate than one running sum
static double dot (const double* x, const double* y, size_t n) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        s0 += x[i] * y[i];
        s1 += x[i + 1] * y[i + 1];
        s2 += x[i + 2] * y[i + 2];
        s3 += x[i + 3] * y[i + 3];
    }
    for (; i < n; i++)
        s0 += x[i] * y[i];
    return (s0 + s1) + (s2 + s3);
}

// apply the plane rotation [x y] <- [x y] * [[c, s], [-s, c]]
static void rotate (double* x, double* y, size_t n, double c, double s) {
    for (size_t i = 0; i < n; i++) {
        double xi = x[i];
        double yi = y[i];
        x[i] = c * xi - s * yi;
        y[i] = s * xi + c * yi;
    }
}

typedef struct SortItem {
    double value;
    int index;
} SortItem;

static int compareDescending (const void* pa, const void* pb) {
    const SortItem* a = (const SortItem*) pa;
    const SortItem* b = (const SortItem*) pb;
    if (a->value > b->value)
        return -1;
    if (a->value < b->value)
        return 1;
    return a->index - b->index;  // stable order for equal values
}

// One-sided Jacobi (Hestenes 1958) on a tall working matrix W (mw x nw, mw >= nw).
// Columns of W are stored contiguously: column j is w->data + j * mw.
// The columns are made mutually orthogonal by plane rotations W <- W * J, while
// the same rotations are accumulated in V (nw x nw, columns contiguous as well).
// At convergence W = U * Sigma, so sigma_j = ||w_j|| and u_j = w_j / sigma_j.
static void jacobiSweeps (Matrix* w, Matrix* v, double* norms, int maxSweeps, JacobiStats* stats) {
    const int nw = w->rows;
    const size_t mw = (size_t) w->cols;
    const size_t nv = (size_t) v->cols;
    // pair (p, q) counts as orthogonal when |<w_p, w_q>| <= tol * ||w_p|| ||w_q||
    const double tol = DBL_EPSILON * (double) mw;

    JacobiStats local = { 0, 0, 0 };

    for (int sweep = 0; sweep < maxSweeps; sweep++) {
        // refresh the squared column norms, they are only updated incrementally below
        double maxNorm = 0.0;
        for (int j = 0; j < nw; j++) {
            const double* wj = w->data + (size_t) j * mw;
            norms[j] = dot (wj, wj, mw);
            if (norms[j] > maxNorm)
                maxNorm = norms[j];
        }
        // columns this small relative to the largest one only carry rounding noise:
        // they correspond to singular values below eps * sigma_max, so leave them alone
        const double negligible = maxNorm * DBL_EPSILON * DBL_EPSILON;

        long rotations = 0;
        for (int p = 0; p < nw - 1; p++) {
            double* wp = w->data + (size_t) p * mw;
            double* vp = v->data + (size_t) p * nv;
            for (int q = p + 1; q < nw; q++) {
                double alpha = norms[p];
                double beta = norms[q];
                if (alpha <= negligible || beta <= negligible)
                    continue;

                double* wq = w->data + (size_t) q * mw;
                double gamma = dot (wp, wq, mw);
                if (fabs (gamma) <= tol * sqrt (alpha) * sqrt (beta))
                    continue;

                // rotation angle that zeroes the (p, q) entry of the 2x2 Gram matrix
                // [[alpha, gamma], [gamma, beta]]; the smaller root of
                // t^2 + 2 zeta t - 1 = 0 keeps |angle| <= pi/4 (Rutishauser)
                double zeta = (beta - alpha) / (2.0 * gamma);
                double t = (zeta >= 0.0 ? 1.0 : -1.0) / (fabs (zeta) + hypot (1.0, zeta));
                double c = 1.0 / sqrt (1.0 + t * t);
                double s = c * t;

                rotate (wp, wq, mw, c, s);
                rotate (vp, v->data + (size_t) q * nv, nv, c, s);

                norms[p] = alpha - t * gamma;
                norms[q] = beta + t * gamma;
                rotations++;
            }
        }

        local.sweeps = sweep + 1;
        local.rotations += rotations;
        if (rotations == 0) {
            local.converged = 1;
            break;
        }
    }

    if (stats != NULL)
        *stats = local;
}

// Householder QR with column pivoting of the tall working matrix W (mw x nw,
// mw >= nw), stored by columns in w (column j is row j of `w`). On return:
//   W P = Q R, column j of W P is column perm[j] of W,
//   R[i][j] (i <= j) is w[j][i],
//   the Householder vector of step j is (vDiag[j], w[j][j + 1 .. mw - 1]),
//   H_j = I - beta[j] v_j v_j^T and Q = H_0 H_1 ... H_{nw-1}.
// Pivoting (Businger, Golub 1965) makes |R_jj| non-increasing, which is what
// makes the subsequent Jacobi iteration on R^T converge quickly.
static void householderQR (Matrix* w, double* vDiag, double* beta, int* perm, double* partial) {
    const int nw = w->rows;
    const size_t mw = (size_t) w->cols;
    double* exact = partial + nw;  // norms at the last exact recomputation

    for (int c = 0; c < nw; c++) {
        const double* y = w->data + (size_t) c * mw;
        partial[c] = dot (y, y, mw);
        exact[c] = partial[c];
        perm[c] = c;
    }

    for (int j = 0; j < nw; j++) {
        const size_t jj = (size_t) j;

        // bring the column with the largest remaining norm to position j
        int pivot = j;
        for (int c = j + 1; c < nw; c++)
            if (partial[c] > partial[pivot])
                pivot = c;
        if (pivot != j) {
            double* pa = w->data + jj * mw;
            double* pb = w->data + (size_t) pivot * mw;
            for (size_t i = 0; i < mw; i++) {
                double t = pa[i];
                pa[i] = pb[i];
                pb[i] = t;
            }
            double t = partial[j];
            partial[j] = partial[pivot];
            partial[pivot] = t;
            t = exact[j];
            exact[j] = exact[pivot];
            exact[pivot] = t;
            int ti = perm[j];
            perm[j] = perm[pivot];
            perm[pivot] = ti;
        }

        double* x = w->data + jj * mw;
        double tail = 0.0;
        for (size_t i = jj + 1; i < mw; i++)
            tail += x[i] * x[i];
        const double norm = sqrt (x[jj] * x[jj] + tail);
        if (tail == 0.0) {
            // already upper triangular in this column: H_j = I
            vDiag[j] = 0.0;
            beta[j] = 0.0;
        } else {
            // reflect x onto -sign (x_j) ||x|| e_j to avoid cancellation in v_j
            const double alpha = x[jj] >= 0.0 ? -norm : norm;
            const double vj = x[jj] - alpha;
            vDiag[j] = vj;
            beta[j] = 2.0 / (vj * vj + tail);
            x[jj] = alpha;

            // apply H_j to the remaining columns
            for (int c = j + 1; c < nw; c++) {
                double* y = w->data + (size_t) c * mw;
                double sum = vj * y[jj];
                for (size_t i = jj + 1; i < mw; i++)
                    sum += x[i] * y[i];
                sum *= beta[j];
                y[jj] -= sum * vj;
                for (size_t i = jj + 1; i < mw; i++)
                    y[i] -= sum * x[i];
            }
        }

        // downdate the remaining column norms; recompute them when cancellation
        // has eaten most of the significant digits (as LAPACK xGEQP3 does)
        for (int c = j + 1; c < nw; c++) {
            const double* y = w->data + (size_t) c * mw;
            partial[c] -= y[jj] * y[jj];
            if (partial[c] <= 1e-6 * exact[c]) {
                double fresh = 0.0;
                for (size_t i = jj + 1; i < mw; i++)
                    fresh += y[i] * y[i];
                partial[c] = fresh;
                exact[c] = fresh;
            }
        }
    }
}

// z <- Q z = H_0 H_1 ... H_{nw-1} z for a vector z of length mw
static void applyQ (const Matrix* w, const double* vDiag, const double* beta, double* z) {
    const size_t mw = (size_t) w->cols;
    for (int j = w->rows - 1; j >= 0; j--) {
        if (beta[j] == 0.0)
            continue;
        const double* x = w->data + (size_t) j * mw;
        const size_t jj = (size_t) j;
        double sum = vDiag[j] * z[jj];
        for (size_t i = jj + 1; i < mw; i++)
            sum += x[i] * z[i];
        sum *= beta[j];
        z[jj] -= sum * vDiag[j];
        for (size_t i = jj + 1; i < mw; i++)
            z[i] -= sum * x[i];
    }
}

// Jacobi SVD with QR preconditioning (Drmac, Veselic 2008; LAPACK xGEJSV):
//   W P = Q R, then one-sided Jacobi on X = R^T (nw x nw): X = U_x S V_x^T, so
//   W = Q X^T P^T = (Q V_x) S (P U_x)^T. Rotating the columns of R^T instead of
//   W is cheaper (columns of length nw instead of mw) and converges in far fewer
//   sweeps, because pivoted QR already orders and concentrates the mass of W.
Svd* svdJacobi (const Matrix* a, int k, int maxSweeps, JacobiStats* stats) {
    const int m = a->rows;
    const int n = a->cols;
    const int minDim = m < n ? m : n;
    if (k <= 0 || k > minDim)
        k = minDim;
    if (maxSweeps <= 0)
        maxSweeps = 40;

    // work with the tall orientation W = A (m >= n) or W = A^T (m < n);
    // the columns of W are the rows of `w`
    const int transposed = m < n;
    const int nw = minDim;
    const int mw = transposed ? n : m;
    Matrix* w = transposed ? copyMatrix (a) : transposeMatrix (a);
    Matrix* x = makeMatrix (nw, nw);   // columns of X = R^T, i.e. the rows of R
    Matrix* vx = makeMatrix (nw, nw);  // columns of V_x
    double* vDiag = (double*) malloc ((size_t) nw * sizeof (double));
    double* beta = (double*) malloc ((size_t) nw * sizeof (double));
    double* norms = (double*) malloc ((size_t) nw * 2 * sizeof (double));
    int* perm = (int*) malloc ((size_t) nw * sizeof (int));
    double* z = (double*) malloc ((size_t) mw * sizeof (double));
    SortItem* order = (SortItem*) malloc ((size_t) nw * sizeof (SortItem));
    Svd* svd = makeSvd (m, n, k);
    Svd* result = NULL;
    if (w == NULL || x == NULL || vx == NULL || vDiag == NULL || beta == NULL ||
        norms == NULL || perm == NULL || z == NULL || order == NULL || svd == NULL)
        goto cleanup;

    householderQR (w, vDiag, beta, perm, norms);
    for (int i = 0; i < nw; i++) {
        for (int j = i; j < nw; j++)
            matrixSet (x, i, j, matrixGet (w, j, i));  // column i of X = row i of R
        matrixSet (vx, i, i, 1.0);
    }

    jacobiSweeps (x, vx, norms, maxSweeps, stats);

    // singular values are the final column norms of X, sorted in decreasing order
    for (int j = 0; j < nw; j++) {
        const double* xj = x->data + (size_t) j * (size_t) nw;
        order[j].value = sqrt (dot (xj, xj, (size_t) nw));
        order[j].index = j;
    }
    qsort (order, (size_t) nw, sizeof (SortItem), compareDescending);

    // numerical rank cut-off: anything below this is rounding noise of size eps * ||A||
    const double sigmaTol = order[0].value * DBL_EPSILON * (double) mw;

    // left singular vectors of W (mw x k) and right ones (nw x k)
    Matrix* wU = transposed ? svd->v : svd->u;
    Matrix* wV = transposed ? svd->u : svd->v;
    for (int r = 0; r < k; r++) {
        const int j = order[r].index;
        const double sigma = order[r].value;
        const double* xj = x->data + (size_t) j * (size_t) nw;
        const double* vj = vx->data + (size_t) j * (size_t) nw;

        // left vector: Q * [v_j; 0]
        memset (z, 0, (size_t) mw * sizeof (double));
        memcpy (z, vj, (size_t) nw * sizeof (double));
        applyQ (w, vDiag, beta, z);
        for (int i = 0; i < mw; i++)
            matrixSet (wU, i, r, z[i]);

        // right vector: P * x_j / sigma_j (zero for a numerically zero singular value)
        if (sigma > sigmaTol && sigma > 0.0) {
            svd->s[r] = sigma;
            for (int i = 0; i < nw; i++)
                matrixSet (wV, perm[i], r, xj[i] / sigma);
        } else {
            svd->s[r] = 0.0;
        }
    }
    result = svd;
    svd = NULL;

cleanup:
    freeMatrix (w);
    freeMatrix (x);
    freeMatrix (vx);
    free (vDiag);
    free (beta);
    free (norms);
    free (perm);
    free (z);
    free (order);
    freeSvd (svd);
    return result;
}

// ---------------------------------------------------------------------------
// randomized SVD

// splitmix64 (Steele, Lea, Flood 2014): tiny, fast and good enough for sketching
static uint64_t nextRandom (uint64_t* state) {
    uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// uniform double in the open interval (0, 1)
static double nextUniform (uint64_t* state) {
    return ((double) (nextRandom (state) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

// fill a matrix with independent N(0, 1) samples (Box-Muller transform)
static void fillGaussian (Matrix* m, uint64_t seed) {
    uint64_t state = seed;
    size_t count = matrixSize (m);
    for (size_t i = 0; i < count; i += 2) {
        double r = sqrt (-2.0 * log (nextUniform (&state)));
        double phi = SVDC_TWO_PI * nextUniform (&state);
        m->data[i] = r * cos (phi);
        if (i + 1 < count)
            m->data[i + 1] = r * sin (phi);
    }
}

Svd* svdRandomized (const Matrix* a, int k, const SvdOptions* opt) {
    const SvdOptions defaults = defaultSvdOptions ();
    if (opt == NULL)
        opt = &defaults;

    const int m = a->rows;
    const int n = a->cols;
    const int minDim = m < n ? m : n;
    if (k <= 0 || k > minDim)
        k = minDim;
    const int p = opt->oversampling > 0 ? opt->oversampling : 0;
    const int q = opt->powerIters > 0 ? opt->powerIters : 0;
    const int l = k + p;

    // the sketch would span the whole column space: nothing to gain, use exact Jacobi
    if (l >= minDim) {
        Svd* full = svdJacobi (a, 0, opt->maxSweeps, NULL);
        if (full != NULL && truncateSvd (full, k) != 0) {
            freeSvd (full);
            return NULL;
        }
        return full;
    }

    Svd* result = NULL;
    Matrix* omega = makeMatrix (n, l);
    Matrix* y = NULL;
    Matrix* z = NULL;
    Matrix* b = NULL;
    Matrix* uFull = NULL;
    Svd* small = NULL;
    if (omega == NULL)
        goto cleanup;

    // stage A: orthonormal basis Q (stored in y) for the range of (A A^T)^q A Omega
    fillGaussian (omega, opt->seed);
    y = multiply (a, omega);
    if (y == NULL || orthonormalizeColumns (y) < 0)
        goto cleanup;
    for (int it = 0; it < q; it++) {
        // re-orthonormalize after every product so that small singular directions
        // are not swamped by rounding errors (subspace iteration, HMT Alg. 4.4)
        z = multiplyTransA (a, y);
        if (z == NULL || orthonormalizeColumns (z) < 0)
            goto cleanup;
        freeMatrix (y);
        y = multiply (a, z);
        freeMatrix (z);
        z = NULL;
        if (y == NULL || orthonormalizeColumns (y) < 0)
            goto cleanup;
    }

    // stage B: B = Q^T A is small (l x n); its exact SVD gives that of Q Q^T A
    b = multiplyTransA (y, a);
    if (b == NULL)
        goto cleanup;
    small = svdJacobi (b, 0, opt->maxSweeps, NULL);
    if (small == NULL)
        goto cleanup;
    uFull = multiply (y, small->u);  // U = Q * U_B
    if (uFull == NULL)
        goto cleanup;

    freeMatrix (small->u);
    small->u = uFull;
    uFull = NULL;
    small->rows = m;
    if (truncateSvd (small, k) != 0)
        goto cleanup;
    result = small;
    small = NULL;

cleanup:
    freeMatrix (omega);
    freeMatrix (y);
    freeMatrix (z);
    freeMatrix (b);
    freeMatrix (uFull);
    freeSvd (small);
    return result;
}

Svd* svdCompute (const Matrix* a, int k, const SvdOptions* opt) {
    if (opt != NULL && opt->method == SVD_JACOBI) {
        Svd* svd = svdJacobi (a, 0, opt->maxSweeps, NULL);
        if (svd != NULL && truncateSvd (svd, k) != 0) {
            freeSvd (svd);
            return NULL;
        }
        return svd;
    }
    return svdRandomized (a, k, opt);
}

int rankForEnergy (const double* s, int count, double totalEnergy, double fraction) {
    if (count <= 0)
        return 0;
    const double target = fraction * totalEnergy;
    double captured = 0.0;
    for (int i = 0; i < count; i++) {
        captured += s[i] * s[i];
        if (captured >= target)
            return i + 1;
    }
    return count;
}

static double capturedEnergy (const Svd* svd, int k) {
    double sum = 0.0;
    for (int i = 0; i < k && i < svd->rank; i++)
        sum += svd->s[i] * svd->s[i];
    return sum;
}

Svd* svdForEnergy (const Matrix* a, double fraction, const SvdOptions* opt) {
    const SvdOptions defaults = defaultSvdOptions ();
    if (opt == NULL)
        opt = &defaults;
    if (fraction > 1.0)
        fraction = 1.0;

    const int minDim = a->rows < a->cols ? a->rows : a->cols;
    const double total = frobeniusNorm2 (a);

    if (opt->method == SVD_JACOBI) {
        Svd* svd = svdJacobi (a, 0, opt->maxSweeps, NULL);
        if (svd == NULL)
            return NULL;
        int k = total > 0.0 ? rankForEnergy (svd->s, svd->rank, total, fraction) : 1;
        if (truncateSvd (svd, k) != 0) {
            freeSvd (svd);
            return NULL;
        }
        return svd;
    }

    // randomized: grow the sketch geometrically until it captures the target energy.
    // randomized singular values never exceed the exact ones, so the chosen k is
    // never smaller than the optimal one (it can be slightly larger)
    int k = minDim < 16 ? minDim : 16;
    for (;;) {
        Svd* svd = svdRandomized (a, k, opt);
        if (svd == NULL)
            return NULL;
        if (total <= 0.0) {
            truncateSvd (svd, 1);
            return svd;
        }
        int reached = capturedEnergy (svd, svd->rank) >= fraction * total;
        if (reached || svd->rank >= minDim) {
            int kk = rankForEnergy (svd->s, svd->rank, total, fraction);
            if (truncateSvd (svd, kk) != 0) {
                freeSvd (svd);
                return NULL;
            }
            return svd;
        }
        freeSvd (svd);
        k = 2 * k < minDim ? 2 * k : minDim;
    }
}

// shrink a rows x oldCols row-major matrix to its first newCols columns, in place
static int dropColumns (Matrix* m, int newCols) {
    const int oldCols = m->cols;
    for (int i = 0; i < m->rows; i++)
        memmove (m->data + (size_t) i * (size_t) newCols,
                 m->data + (size_t) i * (size_t) oldCols,
                 (size_t) newCols * sizeof (double));
    double* shrunk = (double*) realloc (m->data, matrixSize (m) / (size_t) oldCols * (size_t) newCols * sizeof (double));
    if (shrunk != NULL)
        m->data = shrunk;  // a failed shrink keeps the (larger) original block: still valid
    m->cols = newCols;
    return 0;
}

int truncateSvd (Svd* svd, int k) {
    if (k <= 0 || k >= svd->rank)
        return 0;
    dropColumns (svd->u, k);
    dropColumns (svd->v, k);
    double* s = (double*) realloc (svd->s, (size_t) k * sizeof (double));
    if (s != NULL)
        svd->s = s;
    svd->rank = k;
    return 0;
}

Matrix* svdReconstruct (const Svd* svd, int k) {
    if (k <= 0 || k > svd->rank)
        k = svd->rank;
    Matrix* us = makeMatrix (svd->rows, k);
    Matrix* vk = makeMatrix (svd->cols, k);
    Matrix* res = NULL;
    if (us != NULL && vk != NULL) {
        for (int i = 0; i < svd->rows; i++)
            for (int r = 0; r < k; r++)
                matrixSet (us, i, r, matrixGet (svd->u, i, r) * svd->s[r]);
        for (int j = 0; j < svd->cols; j++)
            for (int r = 0; r < k; r++)
                matrixSet (vk, j, r, matrixGet (svd->v, j, r));
        res = multiplyTransB (us, vk);  // (U_k S_k) * V_k^T
    }
    freeMatrix (us);
    freeMatrix (vk);
    return res;
}
