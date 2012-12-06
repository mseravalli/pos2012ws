/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"

#include "initialization.h"
#include "test_functions.h"
#include "util_read_files.h"
#include "util_write_files.h"

#define RECV_ELEM  5
#define SEND_ELEM  10
#define INNER_ELEM 15

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int num_elems, double *cgup_local) {
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
    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf,
                                points_count, points, elems);
    vtk_append_double(file_vtk_out, "DISTRIBUTION", nintci, nintcf, distr);

    free(distr);

    return result;
}

int test_communication(char* file_in, char* file_vtk_out, int* local_global_index,
                       int* num_elems, int neighbors_count,
                       int* send_count, int** send_list,
                       int* recv_count, int** recv_list) {
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

    int* commlist = (int*) malloc(elem_count * sizeof(int));
    for (int i = 0; i < elem_count; ++i) {
        commlist[i] = 0;
    }

    for (int i = 0; i < *num_elems; ++i) {
        commlist[local_global_index[i]] = INNER_ELEM;
    }

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < send_count[i]; ++j) {
            commlist[local_global_index[send_list[i][j]]] = SEND_ELEM;
        }
        for (int j = 0; j < recv_count[i]; ++j) {
            commlist[recv_list[i][j]] = RECV_ELEM;
        }
    }

    // write the result to a vtk file
    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf,
                                points_count, points, elems);
    vtk_append_integer(file_vtk_out, "COMMUNICATION", nintci, nintcf, commlist);

    free(commlist);

    return result;
}



