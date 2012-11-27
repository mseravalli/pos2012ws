/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "initialization.h"
#include "test_functions.h"
#include "util_read_files.h"
#include "util_write_files.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index, 
                      int num_elems, double *cgup_local) {
    // TODO: return proper values for the correct execution
    int result = 0;
    
    int nintci, nintcf;
    int nextci, nextcf;
    int **lcc;
    double *bs, *be, *bn, *bw, *bl, *bh;
    double *bp;
    double *su;
    int points_count;
    int** points;
    int* elems;

    // initialize the elements
    result = read_binary_geo(file_in, 
                             &nintci, &nintcf, &nextci, &nextcf,
                             &lcc,
                             &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, 
                             &points_count, &points, &elems);
    
    if ( result != 0 ) return result;

    // select the correct content to instert into distr
    int elem_count = nintcf -nintci + 1; 

    double* distr = (double*) calloc(elem_count, sizeof(double));

    for (int i = 0; i < elem_count; ++i) {
        distr[i] = 0.0;
    }

    for (int i = 0; i < num_elems; ++i) {
        distr[local_global_index[i]] = cgup_local[i];
    }  

    // write the result to a vtk file
    char file_out[100];
    sprintf(file_out, "partition.vtk");
    vtk_write_unstr_grid_header(file_in, file_out, nintci, nintcf, points_count, 
                                points, elems);
    vtk_append_double(file_out, "DISTRIBUTION", nintci, nintcf, distr);

    free(distr);

    return result;
}

int test_communication(char* file_in, char* file_vtk_out, int* local_global_index,
                       int* num_elems, int neighbors_count, 
                       int* send_count, int** send_list, 
                       int* recv_count, int** recv_list) {
    // TODO: return the appropriate values
    int result = 0;

    int* commlist = (int*) malloc((*num_elems) * sizeof(int));
    
    for (int i = 0; i < (*num_elems); ++i) {
        commlist[i] = 0;
    }

    return result;
}



