#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <memory.h>

#include <lapacke.h>
#include <lapack.h>
#include <cblas.h>

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
            double value = data [m*i + j];
            if (value < 0) {
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
    int size = fmin (row_b, column_b);
	double* B = (double*) malloc (sizeof(double)*size*size);
    int k = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (j == i) { 
				B[i*size+j] = A[k];
                k++;
            }
            else
                B[i*size+j] = 0;
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

double* getMult (double* U, double* S, double* V, int row, int column, int k) {
    double* intermediate = (double*) malloc (sizeof(double)*row*column);
    double* res = (double*) malloc (sizeof(double)*row*column);

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            res [i * column + j]  = 0;
            for (int r = 0; r < k; r++) { 
                intermediate [i*column + j] += U [i * column + r] * S [r * row + j];
            }
        }
    }

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < column; j++) {
            res [i * column + j]  = 0;
            for (int r = 0; r < k; r++) { 
                res [i*column + j] += intermediate [i * column + r] * V [r * row + j];
            }
        }
    }
    
    return res;
}

double* shrinker (double* A, int rows, int cols, int rowsNew, int colsNew) {
    double* res = (double*) malloc (sizeof(double) * rows * cols);
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            if (i<rowsNew && j<colsNew) 
                res [cols*i + j] = A [cols*i + j];
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
	lapack_int i = LAPACKE_dgesvd(LAPACK_ROW_MAJOR,'A', 'A', M, N, A, M, S, U, M, V, N, superb);
    printf ("\t\tLAPACKE returned: %d\n", i);
 	printf ("\t\tLAPACKE function complete.\n");

 	printf ("\t\tAttempting to make S diagonal...\n");
	double* S_ = diag(S, row, column);
 	printf ("\t\tDiagonalization complete.\n");
    return S_;
}

Matrix* svdCompression (Matrix* A, int k) {
    int column = A->width;
    int row = A->height;

    printf ("\tAttempting to allocate matrices for SVD...\n");
    double* matrix = toDouble (A->data, row, column);
    double* U = (double*) malloc (sizeof(double)*row*row);
    double* S = (double*) malloc (sizeof(double)*row*column);
    double* V = (double*) malloc (sizeof(double)*column*column);
    printf ("\tAllocation complete.\n");

    printf ("\tAttempting to apply SVD to the input matrix...\n");
    double* S_ = shrinker(svd (matrix, U, S, V, row, column), row, column, k, k);
    printf ("\tApplying SVD to the input matrix complete.\n");

    //double* U_ = shrinker(U, row, row, row, k);
    //double* V_ = shrinker(V, column, column, k, row);
    
    printf ("\tAttempting to multiply singular matrices...\n");
    double* resArrDouble = getMult (U, S_, V, row, column, k);
    unsigned char* resArrChar = toChar (resArrDouble, row, column);
    printf ("\tMultiplication complete.\n");

    Matrix* res = makeMatrix (column, row, 1);
    dataToMatrix (res, resArrChar, row, column, 1);

    return res;
}