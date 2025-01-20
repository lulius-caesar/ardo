/* 
 * CLASSICAL LAMINATE THEORY 
 * clt.h
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.c"
#include "rw.c"
#include "clt.h"


void 
assemble_Z(double *Z, Laminate *lam, bool debugZ) 
{
        double total_thk = 0.0;
        double h;
        for(i = 0; i < lam->num_layers; i++) h += lam->thk[i];

        for (i = 0; i < lam->num_layers; ++i) {
                total_thk += lam->thk[i];
        }

        Z[0] = -total_thk / 2.0;
        for (i = 1; i <= lam->num_layers; ++i) {
                Z[i] = Z[i-1] + lam->thk[i-1];
        }

        if (debugZ) {
                printf("h = %0.4e mm\n", h);
                printf("Layer Interface Positions (Z):\n");
                for (i = 0; i <= lam->num_layers; ++i) {
                        printf("Z[%d] = %.4e mm\n", i, Z[i]);
                }
        }
}

void 
calc_thermal_forces(double *Nt, MaterialProperties *mat_list, 
                    Laminate *lam, double *Z, int *fail_list, double dT) 
{
        double a[3], a_LCS[3];
        double Q[3][3], Qbar[3][3], T[3][3], Ti[3][3], temp[3][3];

        for (i = 0; i < lam->num_layers; ++i) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];
                a[0] = mat_prop->a1; a[1] = mat_prop->a2; a[2] = 0.0;

                assemble_T(T, lam->ang[i]);

                /* Transform to Laminate Coordinate System */
                for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                                a_LCS[j] += T[k][j] * a[k];
                        }
                }

                /* Assemble Q matrix */
                assemble_Q(Q, mat_prop);

                /* Calculate Qbar */
                transpose(T,Ti);
                matmul(Ti,Q,temp);
                matmul(Q,T,Qbar);

                /* Calculate force resultants */
                for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                                double force = Qbar[j][k] * a_LCS[k] * dT;
                                Nt[j] += force * (Z[i+1] - Z[i]);
                                Nt[j+3] += force * 0.5 * 
                                        (Z[i+1]*Z[i+1] - Z[i]*Z[i]);
                        }
                }
        }
}

void 
calc_moisture_forces(double *Nm, MaterialProperties *mat_list, 
                     Laminate *lam, double *Z, int *fail_list, double dM) 
{
        double b[3], b_LCS[3];
        double Q[3][3], Qbar[3][3], T[3][3], Ti[3][3], temp[3][3];

        for (i = 0; i < lam->num_layers; ++i) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];
                b[0] = mat_prop->b1; b[1] = mat_prop->b2; b[2] = 0.0;

                /* Apply transformed matrix */
                assemble_T(T, lam->ang[i]); 

                /* Transform to Laminate Coordinate System */
                for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                                b_LCS[j] += T[k][j] * b[k];
                        }
                }

                /* Assemble Q matrix */
                assemble_Q(Q, mat_prop);

                /* Calculate Qbar */
                transpose(T,Ti);
                matmul(Ti,Q,temp);
                matmul(Q,T,Qbar);

                /* Calculate force resultants */
                for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                                double force = Qbar[j][k] * b_LCS[k] * dM;
                                Nm[j] += force * (Z[i+1] - Z[i]);
                                Nm[j+3] += force * 0.5 * 
                                        (Z[i+1]*Z[i+1] - Z[i]*Z[i]);
                        }
                }
        }
}

void 
assemble_Q(double Q[3][3], MaterialProperties *mat_prop) 
{

        if (mat_prop->isotropic) {
                double E, nu, E_div;
                E = mat_prop->E1;
                nu = mat_prop->nu12;

                /*if (fail_type == 1 || fail_type == 2 || fail_type == 3) {*/
                /*        double df = 0.001; */
                /*        E *= df;*/
                /*        nu *= df;*/
                /*}*/
                /**/
                E_div = E / (1 - nu * nu);

                Q[0][0] = Q[1][1] = E_div;
                Q[0][1] = Q[1][0] = nu * E_div;
                Q[2][2] = E_div * (1 - nu) / 2.0;  // Since G = E / (2(1+nu)) for isotropic
                Q[0][2] = Q[1][2] = Q[2][0] = Q[2][1] = 0.0;
        } else {
                double E1 = mat_prop->E1, 
                E2 = mat_prop->E2, 
                nu12 = mat_prop->nu12, 
                G12 = mat_prop->G12;
                double df = 0.001;  /* Degradation factor */

                /* fiber or shear failure */
                /*if (fail_type == 1 || fail_type == 2) {  */
                /*        E1 *= df; */
                /*        E2 *= df; */
                /*        nu12 *= df; */
                /*        G12 *= df;*/
                /*} else if (fail_type == 3) {  */
                /* matrix failure */
                /*        E2 *= df; */
                /*        nu12 *= df; */
                /*        G12 *= df;*/
                /*}*/

                double n21 = nu12 * E2 / E1;

                Q[0][0] = E1 / (1 - nu12 * n21);
                Q[0][1] = Q[1][0] = nu12 * E1 * E2 / (E1 - nu12 * nu12 * E2);
                Q[1][1] = E1 * E2 / (E1 - nu12 * nu12 * E2);
                Q[2][2] = G12;
                Q[0][2] = Q[1][2] = Q[2][0] = Q[2][1] = 0.0;
        }
}

void 
assemble_T(double T[3][3], double angle) 
{
        double rad = angle * DEGTORAD;

        double cos = Cos(rad), sin = Sin(rad);
        double cs = cos * sin, cc = cos * cos, ss = sin * sin;

        T[0][0] = cc;    T[0][1] = ss;    T[0][2] = cs;
        T[1][0] = ss;    T[1][1] = cc;    T[1][2] = -cs;
        T[2][0] = -2*cs; T[2][1] = 2*cs;  T[2][2] = cc - ss;
}

void 
assemble_ABD(double ABD[6][6], MaterialProperties *mat_list, 
             Laminate *lam, double *Z, bool debugABD)
{
        double A[3][3] = {0}, B[3][3] = {0}, D[3][3] = {0};
        double Q[3][3], T[3][3], Ti[3][3];
        double temp[3][3], Qbar[3][3];

        for (i = 0; i < lam->num_layers; ++i) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];

                assemble_Q(Q, mat_prop);
                assemble_T(T, lam->ang[i]);

                /* Transpose T to get Ti */
                transpose(T,Ti);

                /* Qbar = T^T * Q * T */
                matmul(Ti, Q, temp);
                matmul(temp, T, Qbar);

                /* B and D matrices */
                //double h = Z[i+1] - Z[i]; /* thickness of the layer */
                for (j = 0; j < 3; ++j) {
                        for (k = 0; k < 3; ++k) {
                                A[j][k] += Qbar[j][k] * (Z[i+1] - Z[i]);
                                B[j][k] += 0.5 * Qbar[j][k] * (Z[i+1] * Z[i+1] - Z[i] * Z[i]);
                                D[j][k] += CONST13 * Qbar[j][k] * (Z[i+1] * Z[i+1] * Z[i+1] - Z[i] * Z[i] * Z[i]);
                        }
                }
        }

        /* Combine into ABD matrix */
        for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                        ABD[i][j] = A[i][j];           /* A matrix */
                        ABD[i][j+3] = ABD[i+3][j] = B[i][j];  /* B matrix */
                        ABD[i+3][j+3] = D[i][j];       /* D matrix */
                }
        }
        if (debugABD) {
                dispmat(&ABD[0][0], 6, 6, 0);
        }
}

/* Tsai-Wu criterion */
void 
fs_tsaiwu_2D(MaterialProperties *mat_prop, 
             double sig1, double sig2, double tau, 
             FSResult *result) 
{
        double f11 = 1 / (mat_prop->Xt * mat_prop->Xc);
        double f22 = 1 / (mat_prop->Yt * mat_prop->Yc);
        double f12 = -1 / (2 * sqrt(mat_prop->Xt * mat_prop->Xc * mat_prop->Yt * mat_prop->Yc));
        double f66 = 1 / (mat_prop->S12 * mat_prop->S12);
        double f1 = 1/mat_prop->Xt - 1/mat_prop->Xc;
        double f2 = 1/mat_prop->Yt - 1/mat_prop->Yc;

        double a = f11*sig1*sig1 + f22*sig2*sig2 + f66*tau*tau + 2*f12*sig1*sig2;
        double b = f1*sig1 + f2*sig2;

        result->fs = (-b + sqrt(b*b + 4*a)) / (2*a);

        double H1 = fabs(f1*sig1 + f11*sig1*sig1);
        double H2 = fabs(f2*sig2 + f22*sig2*sig2);
        double H6 = fabs(f66*tau*tau);

        if (H1 >= H2 && H1 >= H6) {
                strcpy(result->mode, "fiber");
        } else if (H2 >= H1 && H2 >= H6) {
                strcpy(result->mode, "matrix");
        } else {
                strcpy(result->mode, "shear");
        }
}

void 
tsaiwu_2D(MaterialProperties *mat_list, Laminate *lam, double **stress_inf,
          double **stress_sup, FSResult *fs_inf, FSResult *fs_sup) 
{
        for (i = 0; i < lam->num_layers; i++) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];

                /* Superior Stresses */
                fs_tsaiwu_2D(mat_prop, stress_sup[0][i], stress_sup[1][i], stress_sup[2][i], &fs_sup[i]);

                /* Inferior Stresses */
                fs_tsaiwu_2D(mat_prop, stress_inf[0][i], stress_inf[1][i], stress_inf[2][i], &fs_inf[i]);
        }
}

/* Maximum Stress criterion */
void 
fs_maxstress_2D(MaterialProperties *mat_prop, double sig1, double sig2, 
                double tau, FSResult *result) 
{
        double f_1 = sig1 > 0 ? sig1 / mat_prop->Xt : -sig1 / mat_prop->Xc;
        double f_2 = sig2 > 0 ? sig2 / mat_prop->Yt : -sig2 / mat_prop->Yc;
        double f_s = fabs(tau) / mat_prop->S12;

        double f_max = fmax(fmax(f_1, f_2), f_s);

        result->fs = 1 / f_max;

        if (f_max == f_1) {
                strcpy(result->mode, "fiber");
        } else if (f_max == f_2) {
                strcpy(result->mode, "matrix");
        } else {
                strcpy(result->mode, "shear");
        }
}

void 
maxstress_2D(MaterialProperties *mat_list, Laminate *lam, 
             double **stress_inf, double **stress_sup, 
             FSResult *fs_inf, FSResult *fs_sup) 
{
        for (i = 0; i < lam->num_layers; i++) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];

                /* Superior Stresses */
                fs_maxstress_2D(mat_prop, stress_sup[0][i], stress_sup[1][i],
                                stress_sup[2][i], &fs_sup[i]);

                /* Inferior Stresses */
                fs_maxstress_2D(mat_prop, stress_inf[0][i], stress_inf[1][i],
                                stress_inf[2][i], &fs_inf[i]);
        }
}

/* Maximum Strain criterion is not directly translated due to lack of strain input in the given code. */

/* Hashin criterion */
void 
fs_hashin_2D(MaterialProperties *mat_prop, 
             double sig1, double sig2, double tau, 
             FSResult *result) 
{
        double f_1 = sig1 >= 0 ? sig1 / mat_prop->Xt : -sig1 / mat_prop->Xc;
        double f_2 = sig2 >= 0 ? sqrt((sig2/mat_prop->Yt)*(sig2/mat_prop->Yt) 
                                      + (tau/mat_prop->S12)*(tau/mat_prop->S12)) :
                sqrt((sig2/mat_prop->Yc)*(sig2/mat_prop->Yc) + 
                     (tau/mat_prop->S12)*(tau/mat_prop->S12));

        double f_max = fmax(f_1, f_2);

        result->fs = 1 / f_max;

        if (f_max == f_1) {
                strcpy(result->mode, "fiber");
        } else {
                strcpy(result->mode, "matrix");
        }
}

void 
hashin_2D(MaterialProperties *mat_list, Laminate *lam, 
          double **stress_inf, double **stress_sup, 
          FSResult *fs_inf, FSResult *fs_sup) 
{
        for (i = 0; i < lam->num_layers; i++) {
                MaterialProperties *mat_prop = &mat_list[lam->mat_id[i]];

                /* Superior Stresses */
                fs_hashin_2D(mat_prop, stress_sup[0][i], stress_sup[1][i], stress_sup[2][i], &fs_sup[i]);

                /* Inferior Stresses */
                fs_hashin_2D(mat_prop, stress_inf[0][i], stress_inf[1][i], stress_inf[2][i], &fs_inf[i]);
        }
}

double E_x(double ABD[6][6], double h) {
        double matrix1[6][6], matrix2[5][5];
        double det1, det2, Ex;

        // Copy relevant parts of ABD to matrix1 for the first determinant
        for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                        matrix1[i][j] = ABD[i][j]; 
                }
        }

        // Perform LU decomposition for matrix1
        if (LUDecompose(&matrix1[0][0], 6, 6) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return -1;
        }

        // Calculate determinant for matrix1
        det1 = LUDeterminant(&matrix1[0][0], 6, 6);

        matrix2[0][0] = ABD[1][1];
        matrix2[0][1] = ABD[1][2];
        matrix2[0][2] = ABD[1][3];
        matrix2[0][3] = ABD[1][4];
        matrix2[0][4] = ABD[1][5];
        matrix2[1][0] = ABD[2][1];
        matrix2[1][1] = ABD[2][2];
        matrix2[1][2] = ABD[2][3];
        matrix2[1][3] = ABD[2][4];
        matrix2[1][4] = ABD[2][5];
        matrix2[2][0] = ABD[3][1];
        matrix2[2][1] = ABD[3][2];
        matrix2[2][2] = ABD[3][3];
        matrix2[2][3] = ABD[3][4];
        matrix2[2][4] = ABD[3][5];
        matrix2[3][0] = ABD[4][1];
        matrix2[3][1] = ABD[4][2];
        matrix2[3][2] = ABD[4][3];
        matrix2[3][3] = ABD[4][4];
        matrix2[3][4] = ABD[4][5];
        matrix2[4][0] = ABD[5][1];
        matrix2[4][1] = ABD[5][2];
        matrix2[4][2] = ABD[5][3];
        matrix2[4][3] = ABD[5][4];
        matrix2[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix2
        if (LUDecompose(&matrix2[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix2.\n");
                return -1;
        }

        // Calculate determinant for matrix2
        det2 = LUDeterminant(&matrix2[0][0], 5, 5);

        // Check for division by zero
        if (det2 == 0) {
                fprintf(stderr, "Error: Division by zero in E_x calculation, determinant of the submatrix is zero.\n");
                return -1; 
        }

        // Compute E_x
        Ex = (det1 / det2) / h;
        return Ex;
}

void
E_y(double ABD[6][6], double h, double Ey) 
{
        double matrix1[6][6], matrix2[5][5];
        double det1, det2;

        // Copy relevant parts of ABD to matrix1 for the first determinant
        for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                        matrix1[i][j] = ABD[i][j]; 
                }
        }

        // Perform LU decomposition for matrix1
        if (LUDecompose(&matrix1[0][0], 6, 6) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return;
        }

        // Calculate determinant for matrix1
        det1 = LUDeterminant(&matrix1[0][0], 6, 6);

        matrix2[0][0] = ABD[0][0];
        matrix2[0][1] = ABD[0][2];
        matrix2[0][2] = ABD[0][3];
        matrix2[0][3] = ABD[0][4];
        matrix2[0][4] = ABD[0][5];
        matrix2[1][0] = ABD[2][0];
        matrix2[1][1] = ABD[2][2];
        matrix2[1][2] = ABD[2][3];
        matrix2[1][3] = ABD[2][4];
        matrix2[1][4] = ABD[2][5];
        matrix2[2][0] = ABD[3][0];
        matrix2[2][1] = ABD[3][2];
        matrix2[2][2] = ABD[3][3];
        matrix2[2][3] = ABD[3][4];
        matrix2[2][4] = ABD[3][5];
        matrix2[3][0] = ABD[4][0];
        matrix2[3][1] = ABD[4][2];
        matrix2[3][2] = ABD[4][3];
        matrix2[3][3] = ABD[4][4];
        matrix2[3][4] = ABD[4][5];
        matrix2[4][0] = ABD[5][0];
        matrix2[4][1] = ABD[5][2];
        matrix2[4][2] = ABD[5][3];
        matrix2[4][3] = ABD[5][4];
        matrix2[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix2
        if (LUDecompose(&matrix2[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix2.\n");
                return;
        }

        // Calculate determinant for matrix2
        det2 = LUDeterminant(&matrix2[0][0], 5, 5);

        // Check for division by zero
        if (det2 == 0) {
                fprintf(stderr, "Error: Division by zero in E_x calculation, determinant of the submatrix is zero.\n");
                return; 
        }

        // Compute E_x
        Ey = (det1 / det2) / h;
}

double 
G_xy(double ABD[6][6], double h) 
{
        double matrix1[6][6], matrix2[5][5];
        double det1, det2, Gxy;

        // Copy relevant parts of ABD to matrix1 for the first determinant
        for (i = 0; i < 6; i++) {
                for (j = 0; j < 6; j++) {
                        matrix1[i][j] = ABD[i][j]; 
                }
        }

        // Perform LU decomposition for matrix1
        if (LUDecompose(&matrix1[0][0], 6, 6) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return -1;
        }

        // Calculate determinant for matrix1
        det1 = LUDeterminant(&matrix1[0][0], 6, 6);

        matrix2[0][0] = ABD[0][0];
        matrix2[0][1] = ABD[0][1];
        matrix2[0][2] = ABD[0][3];
        matrix2[0][3] = ABD[0][4];
        matrix2[0][4] = ABD[0][5];

        matrix2[1][0] = ABD[1][0];
        matrix2[1][1] = ABD[1][1];
        matrix2[1][2] = ABD[1][3];
        matrix2[1][3] = ABD[1][4];
        matrix2[1][4] = ABD[1][5];

        matrix2[2][0] = ABD[3][0];
        matrix2[2][1] = ABD[3][1];
        matrix2[2][2] = ABD[3][3];
        matrix2[2][3] = ABD[3][4];
        matrix2[2][4] = ABD[3][5];

        matrix2[3][0] = ABD[4][0];
        matrix2[3][1] = ABD[4][1];
        matrix2[3][2] = ABD[4][3];
        matrix2[3][3] = ABD[4][4];
        matrix2[3][4] = ABD[4][5];

        matrix2[4][0] = ABD[5][0];
        matrix2[4][1] = ABD[5][1];
        matrix2[4][2] = ABD[5][3];
        matrix2[4][3] = ABD[5][4];
        matrix2[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix2
        if (LUDecompose(&matrix2[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix2.\n");
                return -1;
        }

        // Calculate determinant for matrix2
        det2 = LUDeterminant(&matrix2[0][0], 5, 5);

        // Check for division by zero
        if (det2 == 0) {
                fprintf(stderr, "Error: Division by zero in E_x calculation, determinant of the submatrix is zero.\n");
                return -1; 
        }

        // Compute E_x
        Gxy = (det1 / det2) / h;
        return Gxy;
}

double 
nu_xy(double ABD[6][6], double h) 
{
        double matrix1[5][5], matrix2[5][5];
        double det1, det2, nxy;

        matrix1[0][0] = ABD[0][1];
        matrix1[0][1] = ABD[1][2];
        matrix1[0][2] = ABD[1][3];
        matrix1[0][3] = ABD[1][4];
        matrix1[0][4] = ABD[1][5];
        matrix1[1][0] = ABD[0][2];
        matrix1[1][1] = ABD[2][2];
        matrix1[1][2] = ABD[2][3];
        matrix1[1][3] = ABD[2][4];
        matrix1[1][4] = ABD[2][5];
        matrix1[2][0] = ABD[0][3];
        matrix1[2][1] = ABD[3][2];
        matrix1[2][2] = ABD[3][3];
        matrix1[2][3] = ABD[3][4];
        matrix1[2][4] = ABD[3][5];
        matrix1[3][0] = ABD[0][4];
        matrix1[3][1] = ABD[4][2];
        matrix1[3][2] = ABD[4][3];
        matrix1[3][3] = ABD[4][4];
        matrix1[3][4] = ABD[4][5];
        matrix1[4][0] = ABD[0][5];
        matrix1[4][1] = ABD[5][2];
        matrix1[4][2] = ABD[5][3];
        matrix1[4][3] = ABD[5][4];
        matrix1[4][4] = ABD[5][5];

        /* Perform LU decomposition for matrix1 */
        if (LUDecompose(&matrix1[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return -1;
        }

        /* Calculate determinant for matrix1 */
        det1 = LUDeterminant(&matrix1[0][0], 5, 5);

        matrix2[0][0] = ABD[1][1];
        matrix2[0][1] = ABD[1][2];
        matrix2[0][2] = ABD[1][3];
        matrix2[0][3] = ABD[1][4];
        matrix2[0][4] = ABD[1][5];
        matrix2[1][0] = ABD[2][1];
        matrix2[1][1] = ABD[2][2];
        matrix2[1][2] = ABD[2][3];
        matrix2[1][3] = ABD[2][4];
        matrix2[1][4] = ABD[2][5];
        matrix2[2][0] = ABD[3][1];
        matrix2[2][1] = ABD[3][2];
        matrix2[2][2] = ABD[3][3];
        matrix2[2][3] = ABD[3][4];
        matrix2[2][4] = ABD[3][5];
        matrix2[3][0] = ABD[4][1];
        matrix2[3][1] = ABD[4][2];
        matrix2[3][2] = ABD[4][3];
        matrix2[3][3] = ABD[4][4];
        matrix2[3][4] = ABD[4][5];
        matrix2[4][0] = ABD[5][1];
        matrix2[4][1] = ABD[5][2];
        matrix2[4][2] = ABD[5][3];
        matrix2[4][3] = ABD[5][4];
        matrix2[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix2
        if (LUDecompose(&matrix2[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for submatrix.\n");
                return -1;
        }

        // Calculate determinant for matrix2
        det2 = LUDeterminant(&matrix2[0][0], 5, 5);

        // Check for division by zero
        if (det2 == 0) {
                fprintf(stderr, "Error: Division by zero in nu_xy calculation, determinant of the submatrix is zero.\n");
                return -1; 
        }

        // Compute E_x
        nxy = det1 / det2;
        return nxy;
}

double 
nu_yx(double ABD[6][6], double h) 
{
        double matrix1[5][5], matrix2[5][5];
        double det1, det2, nyx;

        matrix1[0][0] = ABD[0][1];
        matrix1[0][1] = ABD[0][2];
        matrix1[0][2] = ABD[0][3];
        matrix1[0][3] = ABD[0][4];
        matrix1[0][4] = ABD[0][5];
        matrix1[1][0] = ABD[2][0];
        matrix1[1][1] = ABD[2][2];
        matrix1[1][2] = ABD[2][3];
        matrix1[1][3] = ABD[2][4];
        matrix1[1][4] = ABD[2][5];
        matrix1[2][0] = ABD[3][0];
        matrix1[2][1] = ABD[3][2];
        matrix1[2][2] = ABD[3][3];
        matrix1[2][3] = ABD[3][4];
        matrix1[2][4] = ABD[3][5];
        matrix1[3][0] = ABD[4][0];
        matrix1[3][1] = ABD[4][2];
        matrix1[3][2] = ABD[4][3];
        matrix1[3][3] = ABD[4][4];
        matrix1[3][4] = ABD[4][5];
        matrix1[4][0] = ABD[5][0];
        matrix1[4][1] = ABD[5][2];
        matrix1[4][2] = ABD[5][3];
        matrix1[4][3] = ABD[5][4];
        matrix1[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix1
        if (LUDecompose(&matrix1[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return -1;
        }

        // Calculate determinant for matrix1
        det1 = LUDeterminant(&matrix1[0][0], 5, 5);

        matrix2[0][0] = ABD[0][0];
        matrix2[0][1] = ABD[0][2];
        matrix2[0][2] = ABD[0][3];
        matrix2[0][3] = ABD[0][4];
        matrix2[0][4] = ABD[0][5];
        matrix2[1][0] = ABD[2][0];
        matrix2[1][1] = ABD[2][2];
        matrix2[1][2] = ABD[2][3];
        matrix2[1][3] = ABD[2][4];
        matrix2[1][4] = ABD[2][5];
        matrix2[2][0] = ABD[3][0];
        matrix2[2][1] = ABD[3][2];
        matrix2[2][2] = ABD[3][3];
        matrix2[2][3] = ABD[3][4];
        matrix2[2][4] = ABD[3][5];
        matrix2[3][0] = ABD[4][0];
        matrix2[3][1] = ABD[4][2];
        matrix2[3][2] = ABD[4][3];
        matrix2[3][3] = ABD[4][4];
        matrix2[3][4] = ABD[4][5];
        matrix2[4][0] = ABD[5][0];
        matrix2[4][1] = ABD[5][2];
        matrix2[4][2] = ABD[5][3];
        matrix2[4][3] = ABD[5][4];
        matrix2[4][4] = ABD[5][5];

        // Perform LU decomposition for matrix2
        if (LUDecompose(&matrix2[0][0], 5, 5) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix2.\n");
                return -1;
        }

        // Calculate determinant for matrix2
        det2 = LUDeterminant(&matrix2[0][0], 5, 5);

        // Check for division by zero
        if (det2 == 0) {
                fprintf(stderr, "Error: Division by zero in n_yx calculation, determinant of the submatrix is zero.\n");
                return -1; 
        }

        // Compute E_x
        nyx = det1 / det2;
        return nyx;
}

void 
clt(Laminate *lam, MaterialProperties *materials, 
    int num_materials, double *F, 
    bool print, bool debugABD, bool debugZ) 
{
        /* Allocate and initialize Laminate System strain vectors */
        double LS_strain_sup[3][lam->num_layers];
        double LS_strain_inf[3][lam->num_layers];
        double MS_strain_sup[3][lam->num_layers];
        double MS_strain_inf[3][lam->num_layers];
        double MS_stress_sup[3][lam->num_layers];
        double MS_stress_inf[3][lam->num_layers];
        double T[3][3], Q[3][3];
        double *Z = malloc((lam->num_layers + 1) * sizeof(double));
        assemble_Z(Z, lam, debugZ);

        double ABD[N][N];
        assemble_ABD(ABD, materials, lam, Z, debugABD);

        double h = 0;
        for(i = 0; i < lam->num_layers; i++) h += lam->thk[i];

        if (LUDecompose(&ABD[0][0], 6, 6) < 0) {
                fprintf(stderr, "Error during LU decomposition for matrix1.\n");
                return;
        }
        LUSolve(&ABD[0][0], F, 6, 6);

        double *strain = malloc(3 * sizeof(double));
        double *curvature = malloc(3 * sizeof(double));
        for (i = 0; i < 3; i++) strain[i] = F[i];
        for (i = 0; i < 3; i++) curvature[i] = F[i+3];

        for (i = 0; i < lam->num_layers; i++) {
                for(j = 0; j < 3; j++) {
                        LS_strain_inf[j][i] = strain[j] + curvature[j] * Z[i];
                        LS_strain_sup[j][i] = strain[j] + curvature[j] * Z[i + 1];
                }
        }

        for (i = 0; i < lam->num_layers; i++) {
                MaterialProperties *mat_prop = &materials[lam->mat_id[i]];
                assemble_Q(Q, mat_prop);
                assemble_T(T, lam->ang[i]);
                for (j = 0; j < 3; j++) {
                        MS_strain_sup[j][i] = 0;
                        for (k = 0; k < 3; k++)
                                MS_strain_sup[j][i] += T[j][k] * LS_strain_sup[k][i];
                }
                for (j = 0; j < 3; j++) {
                        MS_strain_inf[j][i] = 0;
                        for (k = 0; k < 3; k++)
                                MS_strain_inf[j][i] += T[j][k] * LS_strain_inf[k][i];
                }
                for (j = 0; j < 3; j++) {
                        MS_stress_sup[j][i] = 0;
                        for (k = 0; k < 3; k++)
                                MS_stress_sup[j][i] += Q[j][k] * MS_strain_sup[k][i];
                }
                for (j = 0; j < 3; j++) {
                        MS_stress_inf[j][i] = 0;
                        for (k = 0; k < 3; k++)
                                MS_stress_inf[j][i] += Q[j][k] * MS_strain_inf[k][i];
                }
        }

        if (print) {
                print_equivalent_properties(ABD, h);
                printf("\n");
                printf("Strain vector:\n");
                for (i = 0; i < 6; i++) {
                        printf("%e\n", F[i]);
                }
                print_lcs_mcs(lam->num_layers, 
                              LS_strain_inf, LS_strain_sup,
                              MS_strain_sup, MS_strain_inf, 
                              MS_stress_sup, MS_stress_inf);
        }

        free(Z);
        free(strain);
        free(curvature);
}

void 
print_lcs_mcs(int num_layers, 
              double LS_strain_inf[3][MAX_LAYERS], 
              double LS_strain_sup[3][MAX_LAYERS], 
              double MS_strain_sup[3][MAX_LAYERS], 
              double MS_strain_inf[3][MAX_LAYERS], 
              double MS_stress_sup[3][MAX_LAYERS], 
              double MS_stress_inf[3][MAX_LAYERS]) 
{
        printf("\n");
        printf("+------------------------------------+\n");
        printf("| Laminate Coordinate System Strains |\n");
        printf("+------------------------------------+\n");
        printf("\n");
        printf("Layer | Inferior (Z-)                  | Superior (Z+)\n");
        printf("------|--------------------------------|----------------------------\n");

        for (j = 0; j < num_layers; j++) {
                printf("%4d  |", j + 1);
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&LS_strain_inf[0][0] + i * num_layers + j));
                }
                printf("  |");
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&LS_strain_sup[0][0] + i * num_layers + j));
                }
                printf("\n");
        }

        printf("\n");
        printf("+------------------------------------+\n");
        printf("| Material Coordinate System Strains |\n");
        printf("+------------------------------------+\n");
        printf("\n");
        printf("Layer | Inferior (Z-)                  | Superior (Z+)\n");
        printf("------|--------------------------------|----------------------------\n");

        for (j = 0; j < num_layers; j++) {
                printf("%4d  |", j + 1);
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&MS_strain_inf[0][0] + i * num_layers + j));
                }
                printf("  |");
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&MS_strain_sup[0][0] + i * num_layers + j));
                }
                printf("\n");
        }

        printf("\n");
        printf("+-----------------------------------+\n");
        printf("| Material Coordinate System Stress |\n");
        printf("+-----------------------------------+\n");
        printf("\n");
        printf("Layer | Inferior (Z-)                  | Superior (Z+)\n");
        printf("------|--------------------------------|----------------------------\n");
        for (j = 0; j < num_layers; j++) {
                printf("%4d  |", j + 1);
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&MS_stress_inf[0][0] + i * num_layers + j));
                }
                printf("  |");
                for (i = 0; i < 3; i++) {
                        printf("%10.2e", *(&MS_stress_sup[0][0] + i * num_layers + j));
                }
                printf("\n");
        }
        printf("\n");
}

void 
print_equivalent_properties(double ABD[6][6], double h) 
{
        double Ex = E_x(ABD, h);
        double Ey = 0; E_y(ABD, h, Ey);
        double Gxy = G_xy(ABD, h);
        double nxy = nu_xy(ABD, h);
        double nyx = nu_yx(ABD, h);

        printf("\n");
        printf("+-----------------------+\n");
        printf("| Engineering Constants |\n");
        printf("+-----------------------+\n");
        printf("\n");

        printf("Equivalent elastic modulus Ex:\t%.4e\n", Ex);
        printf("Equivalent elastic modulus Ey:\t%.4e\n", Ey);
        printf("Equivalent shear modulus Gxy:\t%.4e\n", Gxy);
        printf("Equivalent poisson ratio nu_xy:\t%.3f\n", nxy);
        printf("Equivalent poisson ratio nu_yx:\t%.3f\n", nyx);
}
