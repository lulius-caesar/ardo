/* Author: Julio C. Ramirez-Ceballos
 * Licence: ISC
 * Description:
 * 	This code is used to calculate the equivalent properties of a laminate
 * 	with some configuration.
 *
 * 	INPUT:
 * 	- Material properties
 * 	- Laminate
 * 	- Loads
 *
 * 	MAIN FUNCTION: 
 * 	- calc_stressCLT -- Calculates stress and strain vectors according to
 * 	  the Classical Laminate Theory, at the top and bottom of each layer.
 *
 * 	OUTPUT:
 * 	- Laminate Coord. System (LS) strains, top and bottom of layer.
 * 	- Material Coord. System (MS) stresses and strains, top and bottom of
 * 	  layer;
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clt.c"

void write_output
(const char *filename, const char *content) 
{
        FILE *file = fopen(filename, "w");
        if (file == NULL) {
                printf("Error opening file!\n");
                exit(1);
        }
        fprintf(file, "%s", content);
        fclose(file);
}

void print_usage() {
        printf("Usage: ./main [-p | -w]\n");
        printf("  -p : Print output to console\n");
        printf("  -w : Write output to output.res file\n");
}


int main
(int argc, char *argv[])
{
        if (argc != 2 || (strcmp(argv[1], "-p") != 0 && strcmp(argv[1], "-w") != 0)) {
                print_usage();
                return 1;
        }

        int save_to_file = (strcmp(argv[1], "-w") == 0);

        /* Improve function for saving results and debugging */
        /* Option to save or print output */
        char output[10000]; /* Assuming the output won't exceed this size */
        int output_length = 0;

        /* Redirect stdout to our buffer if we're saving to file */
        if (save_to_file) {
                freopen("output.res", "w", stdout);
        }

        int num_materials;
        MaterialProperties *materials = 
                read_material_properties("materials.dat", &num_materials, 1);

        Laminate lam;
        read_laminate_data("laminate.dat", &lam, 1);

        double F[6];
        read_force_vector(F, "loads.dat");

        /* Perform CLT calculation */
        clt(&lam, materials, num_materials, F, 1, 0, 0);

        /* If we're saving to file, we need to capture the output */
        if (save_to_file) {
                fclose(stdout); /* Close the file stream */
                FILE *file = fopen("output.res", "r");
                if (file == NULL) {
                        printf("Error reading file!\n");
                        exit(1);
                }
                while 
                (fgets(output + output_length, 
                       sizeof(output) - output_length, 
                       file) != NULL) 
                {
                        output_length += strlen(output + output_length);
                }
                fclose(file);
                write_output("output.res", output);
        }

        /* Free allocated memory */
        free(lam.thk);
        free(lam.ang);
        free(lam.mat_id);
        free(materials);

        return 0;
}

