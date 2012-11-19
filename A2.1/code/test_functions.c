/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include "stdlib.h"
#include "test_functions.h"
#include "util_read_files.h"
#include "util_write_files.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index, 
                      int num_elems, double *cgup_local) {
    // Return an error if not implemented
    int result = -1;
    
    int nintci, nintcf;
    int nextci, nextcf;
    int **lcc;
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;
    double *su;
    int points_count;
    int** points;
    int* elems;
    result = read_binary_geo(file_in, 
                             &nintci, &nintcf, &nextci, &nextcf,
                             &lcc,
                             &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, 
                             &points_count, &points, &elems);
    
    if ( result != 0 ) return result;

    double* distr = (double*) calloc(num_elems, sizeof(double));

    for (int i = 0; i < num_elems; ++i) {
        distr[i] = 0.0;
    }

    free(distr);

    return result;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int *num_elems,
                       int neighbors_count, int* send_count, int** send_list, int* recv_count,
                       int** recv_list) {
    // Return an error if not implemented
    return -1;
}



