/*
 * rw.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MAX_mat 100 // Assuming a maximum number of mat
#define MAX_LINE_LENGTH 256 // Assuming lines won't exceed this length
#define N 6

/* Structure for material properties */
typedef struct {
        double E1, E2, nu12, G12, 
        a1, a2, b1, b2,
        Xt, Xc, Yt, Yc, S12, S32;
        bool isotropic;
} MaterialProperties;

/* Structure for laminate properties */
typedef struct {
        int *mat_id;
        double *thk;
        double *ang;
        int num_layers;
} Laminate;

void read_laminate_data
(const char *filename, Laminate *lam, bool print_details) 
{
        FILE *file;
        char line[100]; // Assuming each line won't exceed 100 characters
        int layer_count = 0;
        double h = 0;

        // Attempt to open the file
        file = fopen(filename, "r");
        if (file == NULL) {
                fprintf(stderr, "Error opening file %s\n", filename);
                exit(1);
        }

        // Count the number of lines to determine num_layers
        while (fgets(line, sizeof(line), file) != NULL) {
                if (line[0] != '\n' && line[0] != '#') { // Ignore empty lines or comments
                        layer_count++;
                }
        }

        // Allocate memory for the laminate structure
        lam->num_layers = layer_count;
        lam->mat_id = malloc(layer_count * sizeof(int));
        lam->thk = malloc(layer_count * sizeof(double));
        lam->ang = malloc(layer_count * sizeof(double));

        if (lam->mat_id == NULL || lam->thk == NULL || lam->ang == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                exit(1);
        }

        // Reset file pointer to beginning of file
        rewind(file);

        // Read the actual data
        int i = 0;
        while (fgets(line, sizeof(line), file) != NULL) {
                if (line[0] != '\n' && line[0] != '#') { // Skip empty lines or comments
                        // Temporary variables to hold the parsed data from each line
                        int mat_id;
                        double thickness, angle;

                        // Use sscanf to parse the line
                        if (sscanf(line, "%d %lf %lf", 
                                   &mat_id, &thickness, &angle) == 3) {
                                lam->mat_id[i] = mat_id;
                                lam->thk[i] = thickness;
                                lam->ang[i] = angle;
                                i++;
                        } else {
                                fprintf(stderr, "Error reading line: %s\n", line);
                                // Optionally, you could choose to continue or break here
                        }
                }
        }
        for(int i = 0; i < lam->num_layers; i++) h += lam->thk[i];


        // Print the laminate details if requested
        if (print_details) {
                printf("\n");
                printf("+------------------------+\n");
                printf("| Laminate Configuration |\n");
                printf("+------------------------+\n\n");
                printf("Number of Layers: %d\n", lam->num_layers);
                printf("Laminate total thickness = %.4e\n", h);
                printf("Layer | Thickness (mm) | Angle (°) | Material ID\n");
                printf("------|----------------|-----------|------------\n");
                for (int i = 0; i < lam->num_layers; i++) {
                        printf(" %2d   | %13.3f  | %9.2f | %10d\n",
                               i + 1,
                               lam->thk[i],
                               lam->ang[i],
                               lam->mat_id[i]);
                }
                printf("\n");
        }
        fclose(file);
}

MaterialProperties *read_material_properties
(const char *filename, int *n_mat, bool print_details) 
{
        FILE *file;
        char line[MAX_LINE_LENGTH];
        int material_count = 0;
        MaterialProperties *mat = NULL;

        if (print_details){
                printf("\n");
                printf("+---------------------+\n");
                printf("| Material Properties |\n");
                printf("+---------------------+\n");
                printf("\n");
        }
        file = fopen(filename, "r");
        if (file == NULL) {
                fprintf(stderr, "Error opening file %s\n", filename);
                return NULL;
        }

        // First pass: count the number of mat
        while (fgets(line, sizeof(line), file)) {
                if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
                material_count++;
        }

        // Reset file pointer to beginning
        rewind(file);

        // Allocate memory for mat
        *n_mat = material_count;
        mat = (MaterialProperties *)malloc(
                material_count * sizeof(MaterialProperties));
        if (mat == NULL) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                return NULL;
        }

        // Second pass: read the material properties
        int mat_indx = 0;
        while (fgets(line, sizeof(line), file)) {
                if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

                char type;
                // Check if the material type is specified at the beginning of the line
                if (sscanf(line, "%c", &type) != 1) {
                        fprintf(stderr, "Error reading material type from line: %s", line);
                        continue;
                }

                if (type == 'I') { // 'I' for Isotropic
                        // For isotropic mat, we only need E and nu (Poisson's ratio)
                        if (sscanf(line + 1, "%lf %lf", 
                                   &mat[mat_indx].E1, 
                                   &mat[mat_indx].nu12) != 2) {
                                fprintf(stderr, "Error reading isotropic material properties from line: %s", line);
                                continue;
                        }
                        // Set isotropic flag
                        mat[mat_indx].isotropic = true;
                        // Fill in other properties with default or derived values
                        mat[mat_indx].E2 = mat[mat_indx].E1;  // For isotropic mat, E1 = E2
                        mat[mat_indx].G12 = mat[mat_indx].E1 / (2.0 * (1 + mat[mat_indx].nu12));
                        if (print_details) {
                                printf("ISOTROPIC:\n");
                                printf("Material ID %d Properties:\n", mat_indx);
                                printf("  E  (Pa)  : %.2e\n", mat[mat_indx].E1);
                                printf("  nu       : %.3f\n\n", mat[mat_indx].nu12);
                        }
                } else if (type == 'A') { // 'A' for Anisotropic
                        // Read all properties for anisotropic mat
                        if (sscanf(line + 1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                                   &mat[mat_indx].E1, &mat[mat_indx].E2,
                                   &mat[mat_indx].nu12, &mat[mat_indx].G12,
                                   &mat[mat_indx].a1, &mat[mat_indx].a2,
                                   &mat[mat_indx].b1, &mat[mat_indx].b2,
                                   &mat[mat_indx].Xt, &mat[mat_indx].Xc,
                                   &mat[mat_indx].Yt, &mat[mat_indx].Yc,
                                   &mat[mat_indx].S12, &mat[mat_indx].S32) != 14) {
                                fprintf(stderr, "Error reading anisotropic material properties from line: %s", line);
                                continue;
                        }
                        // Set isotropic flag to false
                        mat[mat_indx].isotropic = false;

                        if (print_details) {
                                printf("ORTHOTROPIC:\n");
                                printf("Material ID %d Properties:\n", mat_indx);
                                printf("  E1   (Pa)  : %.2e\n", mat[mat_indx].E1);
                                printf("  E2   (Pa)  : %.2e\n", mat[mat_indx].E2);
                                printf("  nu12       : %.3f\n", mat[mat_indx].nu12);
                                printf("  G12  (Pa)  : %.2e\n", mat[mat_indx].G12);
                                printf("  a1   (1/K) : %.2e\n", mat[mat_indx].a1);
                                printf("  a2   (1/K) : %.2e\n", mat[mat_indx].a2);
                                printf("  b1         : %.2e\n", mat[mat_indx].b1);
                                printf("  b2         : %.2e\n", mat[mat_indx].b2);
                                printf("  Xt   (Pa)  : %.2e\n", mat[mat_indx].Xt);
                                printf("  Xc   (Pa)  : %.2e\n", mat[mat_indx].Xc);
                                printf("  Yt   (Pa)  : %.2e\n", mat[mat_indx].Yt);
                                printf("  Yc   (Pa)  : %.2e\n", mat[mat_indx].Yc);
                                printf("  S12  (Pa)  : %.2e\n", mat[mat_indx].S12);
                                printf("  S23  (Pa)  : %.2e\n\n", mat[mat_indx].S32);
                        }
                } else {
                        fprintf(stderr, "Unknown material type in line: %s", line);
                        continue;
                }
                mat_indx++;
        }
        fclose(file);
        return mat;
}

void read_force_vector(double *F, const char *filename) {
        FILE *file;
        file = fopen(filename, "r");
        if (file == NULL) {
                fprintf(stderr, "Error opening file %s\n", filename);
                exit(1);
        }

        for (int i = 0; i < N; i++) {
                if (fscanf(file, "%lf", &F[i]) != 1) {
                        fprintf(stderr, "Error reading force vector from file %s\n", filename);
                        fclose(file);
                        exit(1);
                }
        }

        fclose(file);
}
