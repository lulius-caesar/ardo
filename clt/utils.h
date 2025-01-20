#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <stdbool.h>

#define PI 3.14159265358979323846
#define DEGTORAD (PI/180)
#define CONST13 (1.0/3.0)
#define MAT_ELEM(mat, row, col) mat[row][col]
#define SMALLEST_PIVOT (1e-5)
#define SIZE 6

/* Function prototypes */
static double factorial(int n);
double Sin(double x);
double Cos(double x);
void matmul(double A[3][3], double B[3][3], double result[3][3]);
void dot_product(double result[3], double A[3][3], double B[3]);
void transpose(double A[3][3], double Ai[3][3]);
void dispmat(double *A, int rows, int cols, bool transpose);
int LUDecompose(double *amat, int n, int numcols);
void LUSolve(double *amat, double *b, int n, int numcols);
double LUDeterminant(double *amat, int n, int numcols);

#endif /* UTILS_H */
