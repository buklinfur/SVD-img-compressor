// helper for tests/numpy_check.py: reads a matrix, writes its SVD.
//
//   svd_tool <jacobi|randomized> <k> <input.bin> <output.bin>
//
// input:  int32 rows, int32 cols, rows*cols float64 (row-major), native byte order
// output: int32 rank, int32 rows, int32 cols, s[rank], U (rows x rank), V (cols x rank)

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"
#include "svd.h"

int main (int argc, char** argv) {
    if (argc != 5) {
        fprintf (stderr, "usage: svd_tool <jacobi|randomized> <k> <input.bin> <output.bin>\n");
        return 2;
    }
    SvdOptions opt = defaultSvdOptions ();
    opt.method = strcmp (argv[1], "jacobi") == 0 ? SVD_JACOBI : SVD_RANDOMIZED;
    int k = atoi (argv[2]);

    FILE* in = fopen (argv[3], "rb");
    if (in == NULL) {
        perror (argv[3]);
        return 1;
    }
    int32_t dims[2];
    if (fread (dims, sizeof (int32_t), 2, in) != 2) {
        fclose (in);
        return 1;
    }
    Matrix* a = makeMatrix (dims[0], dims[1]);
    if (a == NULL || fread (a->data, sizeof (double), matrixSize (a), in) != matrixSize (a)) {
        fclose (in);
        freeMatrix (a);
        return 1;
    }
    fclose (in);

    Svd* svd = svdCompute (a, k, &opt);
    if (svd == NULL) {
        freeMatrix (a);
        return 1;
    }

    FILE* out = fopen (argv[4], "wb");
    if (out == NULL) {
        perror (argv[4]);
        freeSvd (svd);
        freeMatrix (a);
        return 1;
    }
    int32_t head[3] = { svd->rank, svd->rows, svd->cols };
    fwrite (head, sizeof (int32_t), 3, out);
    fwrite (svd->s, sizeof (double), (size_t) svd->rank, out);
    fwrite (svd->u->data, sizeof (double), matrixSize (svd->u), out);
    fwrite (svd->v->data, sizeof (double), matrixSize (svd->v), out);
    fclose (out);

    freeSvd (svd);
    freeMatrix (a);
    return 0;
}
