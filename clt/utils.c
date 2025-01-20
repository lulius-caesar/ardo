#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

/* Factorial for Taylor series */
static double factorial(int n) {
        if (n == 0 || n == 1) return 1;
        return n * factorial(n - 1);
}

/* Sine function using Taylor series */
double Sin
(double x) 
{
        double sum = 0;
        int n;
        x = fmod(x, 2 * PI); /* Reduce x to [0, 2π) */
        for (n = 0; n < 10; ++n) {
                sum += pow(-1, n) * pow(x, 2*n+1) / factorial(2*n+1);
        }
        return sum;
}

/* Cosine function using Taylor series */
double 
Cos(double x) {
        double sum = 0;
        int n;
        x = fmod(x, 2 * PI); /* Reduce x to [0, 2π) */
        for (n = 0; n < 10; ++n) {
                sum += pow(-1, n) * pow(x, 2*n) / factorial(2*n);
        }
        return sum;
}

/* Matrix multiplication */
void 
matmul(double A[3][3], double B[3][3], double result[3][3]) 
{
        for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                        result[i][j] = 0;
                        for (int k = 0; k < 3; k++) {
                                result[i][j] += A[i][k] * B[k][j];
                        }
                }
        }
}

/* Vector-matrix multiplication */
void 
dot_product(double result[3], double A[3][3], double B[3]) 
{
        for(int i = 0; i < 3; i++) {
                result[i] = 0;
                for(int j = 0; j < 3; j++) {
                        result[i] += A[i][j] * B[j];
                }
        }
}

/* Matrix transposition */
void 
transpose(double A[3][3], double Ai[3][3]) 
{
        for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                        Ai[j][k] = A[k][j];
                }
        }
}

/* Display matrix in row or column major order */
void 
dispmat(double *A, int rows, int cols, bool transpose) 
{
        if (transpose) {
                for (int j = 0; j < cols; j++) {
                        printf("%4d  ", j + 1);
                        for (int i = 0; i < rows; i++) {
                                printf("%10.2e", *(A + i * cols + j));
                        }
                        printf("\n");
                }
        } else {
                for (int i = 0; i < rows; i++) {
                        for (int j = 0; j < cols; j++) {
                                printf("%10.2e", *(A + i * cols + j));
                        }
                        printf("\n");
                }
        }
}

/* Static variables for pivot */
static int *pivot = NULL;

/* LU Decomposition */
int LUDecompose(double *amat, int n, int numcols) {
        if (pivot != NULL) free(pivot);
        pivot = malloc(n * sizeof(int));
        if (!pivot) {
                fprintf(stderr, "Error in LUDecompose - malloc\n");
                return -2;
        }

        for (int i = 0; i < n; i++) {
                pivot[i] = i;
        }

        int sign = 1;
        double *dptr1, *dptr2, dtmp1;

        for (int i = 0; i < n-1; i++) {
                dptr1 = amat + i * numcols + i;
                int k = i;
                for (int j = i + 1; j < n; j++) {
                        dptr2 = amat + j * numcols + i;
                        if (fabs(*dptr2) > fabs(*dptr1)) {
                                dptr1 = dptr2;
                                k = j;
                        }
                }

                if (fabs(*dptr1) < SMALLEST_PIVOT) {
                        fprintf(stderr, "Error in LUDecompose - matrix is singular\n");
                        return -1;
                }

                if (k != i) {
                        int temp = pivot[i]; pivot[i] = pivot[k]; pivot[k] = temp;
                        sign *= -1;
                        for (int j = 0; j < n; j++) {
                                dtmp1 = amat[i * numcols + j];
                                amat[i * numcols + j] = amat[k * numcols + j];
                                amat[k * numcols + j] = dtmp1;
                        }
                }

                for (int j = i+1; j < n; j++) {
                        dtmp1 = amat[j * numcols + i] / *dptr1;
                        amat[j * numcols + i] = dtmp1;
                        for (int k = i+1; k < n; k++) {
                                amat[j * numcols + k] -= dtmp1 * amat[i * numcols + k];
                        }
                }
        }
        /*pivot[n-1] = sign;*/
        *(pivot+n-1) = sign;
        return 1;
}

/* Solve system with LU decomposition */
void LUSolve(double *amat, double *b, int n, int numcols) {
        double dtmp1;
        int i, j;

        for (i = 0; i < n; i++) {
                j = pivot[i];
                if (j != i) {
                        dtmp1 = b[i]; b[i] = b[j]; b[j] = dtmp1;
                }
        }

        for (i = 0; i < n; i++) {
                dtmp1 = b[i];
                for (j = 0; j < i; j++) {
                        dtmp1 -= amat[i * numcols + j] * b[j];
                }
                b[i] = dtmp1;
        }

        for (i = n-1; i >= 0; i--) {
                dtmp1 = b[i];
                for (j = n-1; j > i; j--) {
                        dtmp1 -= amat[i * numcols + j] * b[j];
                }
                b[i] = dtmp1 / amat[i * numcols + i];
        }
}

/* Compute determinant from LU decomposition */
double LUDeterminant(double *amat, int n, int numcols) {
        double det = 1.0;
        for (int i = 0; i < n; i++) {
                det *= amat[i * (numcols + 1)];
        }
        return det * pivot[n - 1];
}
