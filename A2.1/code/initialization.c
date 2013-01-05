/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#include "util_read_files.h"
#include "initialization.h"

int dual_partition(char* part_type, 
              int elems_count, int points_count, 
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    idx_t ne = elems_count;
    idx_t nn = points_count;

    idx_t* eptr = malloc((ne + 1) * sizeof(idx_t));
    idx_t* eind = malloc(ne * 8 * sizeof(idx_t));
    idx_t* vwgt = NULL;
    idx_t* vsize = NULL;
    real_t* tpwgts = NULL;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    *epart = malloc(ne * sizeof(idx_t));
    *npart = malloc(nn * sizeof(idx_t));

    // init eprt and eind
    for (int i = 0; i < ne + 1; ++i) {
        eptr[i] = (idx_t) (i * 8);
    }
    for (int i = 0; i < ne * 8; ++i) {
        eind[i] = (idx_t) elems[i];
    }
    
    result = METIS_PartMeshDual(&ne, &nn,
                                eptr, eind,
                                vwgt, vsize,
                                &ncommon, &nparts,
                                tpwgts, options,
                                objval, *epart, *npart);

    if (result == METIS_OK) {
        result = 0;
    }

    free(eptr);
    free(eind);

    return result;
}

int classical_partition(char* part_type, 
                        int elems_count, int points_count, 
                        int* elems,
                        idx_t ncommon, idx_t nparts,
                        idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = 0;
    idx_t ne = elems_count;
    idx_t nn = points_count;

    *epart = malloc(ne * sizeof(idx_t));
    *npart = malloc(nn * sizeof(idx_t));

    for (int i = 0; i < ne; ++i) {
        (*epart)[i] = i % nparts; 
    }

    // TODO: implement the distribution of the points?

    return result;
}

int map_local_global(int elems_count, idx_t* epart, int** local_global,
                     int* local_elems) {
    int result = 0;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *local_elems = 0;
    for (int i = 0; i < elems_count; ++i) {
        if (epart[i] == rank) {
            (*local_elems)++;
        }
    }

    /** 
     * the local_global array will reference the position
     * of the current element within the elems array / 8
     */
    *local_global = malloc((*local_elems) * sizeof(int));
    for (int i = 0, j = 0; i < elems_count; ++i) {
        if (epart[i] == rank) {
            (*local_global)[j] = i;
            ++j;
        }
    }

    return result;
}

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, idx_t** epart, idx_t** npart, idx_t* objval, 
                   int* local_elems ) {
    /********** START INITIALIZATION **********/
    int i = 0;
    // read-in the input file
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                   &*points, &*elems);

    if ( f_status != 0 ) return f_status;

    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ )
        (*cgup)[i] = 1.0 / ((*bp)[i]);

    /** Partition data */

    int part_result = -1;

    idx_t ncommon = 4;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    idx_t nparts = (idx_t) size;

    int elems_count = (*nintcf - *nintci) + 1;
    
    if (strcmp(part_type, "dual") == 0) {
        part_result = dual_partition(part_type, 
                                     elems_count, *points_count, 
                                     *elems,
                                     ncommon, nparts,
                                     objval, epart, npart);
    } else {
        part_result = classical_partition(part_type, 
                                          elems_count, *points_count, 
                                          *elems,
                                          ncommon, nparts,
                                          objval, epart, npart);
    }

    if (part_result != 0) {
        printf("partition failed\n");
        return -1;
    }

    map_local_global(elems_count, *epart, local_global_index, local_elems);

    return 0;
}

