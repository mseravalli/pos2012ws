/**
 * Initialization step - parse the input file, compute data distribution, 
 * initialize LOCAL computational arrays
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

    *epart = malloc(elems_count * sizeof(idx_t));
    *npart = malloc(points_count * sizeof(idx_t));

    int elems_per_part = (int) elems_count / nparts;

    for (int i = 0; i < elems_count; ++i) {
        if ((i / elems_per_part) < nparts) {
            (*epart)[i] = (i / elems_per_part); 
        } else {
            (*epart)[i] = nparts - 1; 
        }
    }

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

/**
 * Initialize the communication list. The communication list will be used
 * subsequently for creating the proper send and receive list
 */
int init_commlist(int local_elems, int* local_global_index, // in, in
                  int elems_count, idx_t* epart, int** lcc, // in, in, in
                  int** commlist, int* neighbors_count,     // out, out
                  int** send_count, int** recv_count) {     // out, out
    int result = 0;
    int my_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // init commlist to 0
    *commlist = (int*) malloc(elems_count * sizeof(int));
    for (int i = 0; i < elems_count; ++i) {
        (*commlist)[i] = 0;
    }

    // keep track of what partitions are my neighbours
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    *send_count = calloc(size, sizeof(int));
    *recv_count = calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        (*send_count)[i] = 0;
        (*recv_count)[i] = 0;
    }

    int e = -1;
    int local_nbr = 0;
    for (int i = 0; i < local_elems; ++i) {
        e = local_global_index[i];
        /**
         * see what cells need to be received and sent
         * if a cell has a neighbour, the neighbour has to be received
         * if a cell has a neighbour, the cell has to be sent
         * exclude the external cells
         */
        for (int j = 0; j < 6; ++j) { 
            int l = lcc[e][j];
            if (l < elems_count) {
                int p = (int) epart[lcc[e][j]]; 
                if (p != my_rank) {
                    (*commlist)[l] = RECV_ELEM;
                    ++((*recv_count)[epart[l]]);

                    (*commlist)[e] = SEND_ELEM;
                    ++((*send_count)[epart[l]]);

                    ++local_nbr;
                }
            }
        }
        // if there are not external neighbours the cell is an internal cell
        if (local_nbr == 0) {
            (*commlist)[e] = INNER_ELEM;
        }
        local_nbr = 0;
    }

    // count neighbours
    (*neighbors_count) = 0;
    for (int i = 0; i < size; ++i) {
        if ((*send_count)[i] > 0) {
            ++(*neighbors_count);
        }
    }
    
    return result;
}

int initialization(char* file_in, char* part_type, 
                   int* nintci, int* nintcf, int* nextci, int* nextcf, 
                   int*** lcc, 
                   double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, 
                   int* points_count, int*** points, 
                   int** elems, 
                   double** var, double** cgup, double** oc, double** cnorm, 
                   int** local_global_index, int** global_local_index,
                   int* neighbors_count, 
                   int** send_count, int*** send_list, 
                   int** recv_count, int*** recv_list, 
                   idx_t** epart, idx_t** npart, idx_t* objval, int* local_elems) {

    /********** START INITIALIZATION **********/
    int i = 0;
    // read-in the input file
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                                   &*points, &*elems);

    if ( f_status != 0 ) return f_status;

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
    
    *global_local_index = calloc(elems_count, sizeof(int));
    for (int i = 0; i < elems_count; ++i) {
        (*global_local_index)[i] = (*epart)[i];    
    }
    map_local_global(elems_count, *epart, local_global_index, local_elems);

    // initialize the arrays
    *var = (double*) calloc(sizeof(double), (*local_elems));
    *cgup = (double*) calloc(sizeof(double), (*local_elems));
    *oc = (double*) calloc(sizeof(double), (*local_elems));
    *cnorm = (double*) calloc(sizeof(double), (*local_elems));

    for (i = 0; i <= 10; i++) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for (i = 0; i < (*local_elems); i++) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }

    for (i = 0; i < (*local_elems); i++) {
        (*cgup)[i] = 1.0 / ((*bp)[(*local_global_index)[i]]);
    }


    /** set up communication */
    int* commlist = NULL;    
    init_commlist(*local_elems, *local_global_index,
                  elems_count, *epart, *lcc,
                  &commlist, neighbors_count,
                  send_count, recv_count);

    free(commlist);

    return 0;
}

