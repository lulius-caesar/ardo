#define MAX_LAYERS 100  /* Assume a maximum number of layers for simplicity */
#define MAX_MODE_LENGTH 10
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))
#define N 6 // Example for a 6x6 matrix, adjust for your needs

static int i, j, k;

/* Custom error for laminate layup issues */
typedef struct {
        char *message;
} LaminateLayupError;

typedef struct {
        double *stress_inf, *stress_sup, *strain_inf, *strain_sup;
} StressStrainResult;

typedef struct {
        double fs;
        char *mode;
} FailureResult;

typedef struct {
        char mode[MAX_MODE_LENGTH];
        double fs;
} FSResult;

/*void calc_stressCLT(MaterialProperties *mat_list, Laminate *lam, */
/*                double F[6], int *fail_list, double dT, double dM);*/
void assemble_Z(double *Z, Laminate *lam, bool debugZ);
void calc_thermal_forces(double *Nt, MaterialProperties *mat_list, 
                Laminate *lam, double *Z, int *fail_list, double dT);
void calc_moisture_forces(double *Nm, MaterialProperties *mat_list,
                Laminate *lam, double *Z, int *fail_list, double dM);
void assemble_Q(double Q[3][3], MaterialProperties *mat_prop);
void assemble_T(double T[3][3], double angle);
void assemble_ABD(double ABD[6][6], MaterialProperties *mat_list, 
                Laminate *lam, double *Z, bool debugABD);
void fs_hashin_2D(MaterialProperties *mat_prop, double sig1, double sig2, 
                double tau, FSResult *result);
void hashin_2D(MaterialProperties *mat_list, Laminate *lam, 
                double **stress_inf, double **stress_sup, 
                FSResult *fs_inf, FSResult *fs_sup);
void fs_tsaiwu_2D(MaterialProperties *mat_prop, 
                double sig1, double sig2, double tau, 
                FSResult *result);
void tsaiwu_2D(MaterialProperties *mat_list, Laminate *lam, 
                double **stress_inf, double **stress_sup, 
                FSResult *fs_inf, FSResult *fs_sup);
void fs_maxstress_2D(MaterialProperties *mat_prop, double sig1, double sig2, 
                double tau, FSResult *result);
void maxstress_2D(MaterialProperties *mat_list, Laminate *lam, 
                double **stress_inf, double **stress_sup, 
                FSResult *fs_inf, FSResult *fs_sup);
double E_x(double ABD[6][6], double h);
void E_y(double ABD[6][6], double h, double Ey);
double G_xy(double ABD[6][6], double h);
double nu_xy(double ABD[6][6], double h);
double nu_yx(double ABD[6][6], double h);
void print_lcs_mcs(int num_layers, 
              double LS_strain_inf[3][MAX_LAYERS], 
              double LS_strain_sup[3][MAX_LAYERS], 
              double MS_strain_sup[3][MAX_LAYERS], 
              double MS_strain_inf[3][MAX_LAYERS], 
              double MS_stress_sup[3][MAX_LAYERS], 
              double MS_stress_inf[3][MAX_LAYERS]);
void print_equivalent_properties(double ABD[6][6], double h);
/* ADD FUNCTION TO COMPUTE THE PRINCIPAL STRESS/STRAIN FOR EACH PLY */
