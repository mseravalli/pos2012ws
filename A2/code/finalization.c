/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "./finalization.h"
#include "./util_write_files.h"

int merge_arrayd(int length_loc, double* array_loc, int* local_global,
                 double** array_glob) {
    int res = 0;

    int size = 0;
    int rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int length_glob = 0;
    MPI_Reduce(&length_loc, &length_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // send all the information to the master process
    MPI_Status status;
    if (rank == 0) {
        *array_glob = (double*) calloc(length_glob, sizeof(double));

        // put the elements of the master in array_glob
        for (int i = 0; i < length_loc; ++i) {
            (*array_glob)[local_global[i]] = array_loc[i];
        }

        for (int i = 1; i < size; ++i) {
            int length_part = 0;
            double* array_part = NULL;
            int* local_global_part = NULL;

            MPI_Recv(&length_part, 1, MPI_INT, i, TAG_LEN, MPI_COMM_WORLD, &status);
            array_part = (double*) calloc(length_part, sizeof(double));
            local_global_part = (int*) calloc(length_part, sizeof(int));

            MPI_Recv(array_part, length_part, MPI_DOUBLE,
                     i, TAG_ARR, MPI_COMM_WORLD, &status);
            MPI_Recv(local_global_part, length_part, MPI_INT,
                     i, TAG_LG, MPI_COMM_WORLD, &status);

            for (int j = 0; j < length_part; ++j) {
                (*array_glob)[local_global_part[j]] = array_part[j];
            }

            free(array_part);
            free(local_global_part);
        }
    } else {
        MPI_Send(&length_loc, 1, MPI_INT, 0, TAG_LEN, MPI_COMM_WORLD);
        MPI_Send(array_loc, length_loc, MPI_DOUBLE, 0, TAG_ARR, MPI_COMM_WORLD);
        MPI_Send(local_global, length_loc, MPI_INT, 0, TAG_LG, MPI_COMM_WORLD);
    }

    return res;
}

void finalization(char* file_in, char* out_prefix, int total_iters,
                  double residual_ratio, int nintci, int nintcf,
                  int points_count, int** points, int* elems, double* var,
                  double* cgup, double* su, int* local_global) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double* var_glob  = NULL;
    double* cgup_glob = NULL;
    double* su_glob   = NULL;
    merge_arrayd(nintcf, var,  local_global, &var_glob);
    merge_arrayd(nintcf, cgup, local_global, &cgup_glob);
    merge_arrayd(nintcf, su,   local_global, &su_glob);

    int el_int_glob = 0;
    MPI_Reduce(&nintcf, &el_int_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        --el_int_glob;
        char file_out[100];
        sprintf(file_out, "%s_summary.out", out_prefix);

        int status = store_simulation_stats(file_in, file_out, nintci,
                                            el_int_glob, var_glob, total_iters,
                                            residual_ratio);

        sprintf(file_out, "%s_data.vtk", out_prefix);
        vtk_write_unstr_grid_header(file_in, file_out,
                                    nintci, el_int_glob,
                                    points_count, points, elems);
        vtk_append_double(file_out, "CGUP", nintci, el_int_glob, cgup_glob);
        vtk_append_double(file_out, "SU",   nintci, el_int_glob, su_glob);
        vtk_append_double(file_out, "VAR",  nintci, el_int_glob, var_glob);

        if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);

        free(var_glob);
        free(cgup_glob);
        free(su_glob);
    }
}

