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
              int el_int_tot, int points_count,
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    idx_t ne = el_int_tot;
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

int nodal_partition(char* part_type,
              int el_int_tot, int points_count,
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    idx_t ne = el_int_tot;
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

    result = METIS_PartMeshNodal(&ne, &nn,
                                 eptr, eind,
                                 vwgt, vsize,
                                 &nparts,
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
                        int el_int_tot, int points_count,
                        int* elems,
                        idx_t ncommon, idx_t nparts,
                        idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = 0;

    *epart = malloc(el_int_tot * sizeof(idx_t));
    *npart = malloc(points_count * sizeof(idx_t));

    int elems_per_part = (int) el_int_tot / nparts;

    for (int i = 0; i < el_int_tot; ++i) {
        if ((i / elems_per_part) < nparts) {
            (*epart)[i] = (i / elems_per_part);
        } else {
            (*epart)[i] = nparts - 1;
        }
    }

    return result;
}

int partition(char* file_in, char* part_type,
              int* nintci, int* nintcf, int* nextci, int* nextcf,
              int*** lcc,
              double** bs, double** be, double** bn, double** bw,
              double** bl, double** bh, double** bp, double** su,
              int* points_count, int*** points,
              int* el_int_tot, int** elems,
              int** global_local_index,
              idx_t** epart, idx_t** npart, idx_t* objval) {
    // read-in the input file
    int f_status = read_binary_geo(file_in,
                                   &*nintci, &*nintcf, &*nextci, &*nextcf,
                                   &*lcc,
                                   &*bs, &*be, &*bn, &*bw,
                                   &*bl, &*bh, &*bp, &*su,
                                   &*points_count, &*points, &*elems);

    if ( f_status != 0 ) return f_status;
    int part_result = -1;

    idx_t ncommon = 4;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    idx_t nparts = (idx_t) size;

    *el_int_tot = (*nintcf - *nintci) + 1;

    if (strcmp(part_type, "dual") == 0) {
        part_result = dual_partition(part_type,
                                     *el_int_tot, *points_count,
                                     *elems,
                                     ncommon, nparts,
                                     objval, epart, npart);
    } else if (strcmp(part_type, "nodal") == 0) {
        part_result = nodal_partition(part_type,
                                      *el_int_tot, *points_count,
                                      *elems,
                                      ncommon, nparts,
                                      objval, epart, npart);
    } else {
        part_result = classical_partition(part_type,
                                          *el_int_tot, *points_count,
                                          *elems,
                                          ncommon, nparts,
                                          objval, epart, npart);
    }

    if (part_result != 0) {
        printf("partition failed\n");
        return -1;
    }

    *global_local_index = calloc(*el_int_tot, sizeof(int));
    for (int i = 0; i < *el_int_tot; ++i) {
        (*global_local_index)[i] = (*epart)[i];
    }

    free(*epart);
    free(*npart);

    return part_result;
}

int map_local_global(int el_int_tot, int* global_local, int** local_global,
                     int* el_int_loc) {
    int result = 0;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    *el_int_loc = 0;
    for (int i = 0; i < el_int_tot; ++i) {
        if (global_local[i] == rank) {
            (*el_int_loc)++;
        }
    }

    /**
     * the local_global array will reference the position
     * of the current element within the elems array / 8
     */
    *local_global = malloc((*el_int_loc) * sizeof(int));
    for (int i = 0, j = 0; i < el_int_tot; ++i) {
        if (global_local[i] == rank) {
            (*local_global)[j] = i;
            ++j;
        }
    }

    return result;
}

/**
 * Initialize the communication lists
 */
int init_commlist(int el_int_loc, int* local_global_index,             // i, i
                  int el_int_tot, int* global_local_index, int** lcc,  // i,i,i
                  int* neighbors_count,                                 // o
                  int** send_count, int*** send_list,                   // o, o
                  int** recv_count, int*** recv_list) {                 // o, o
    int result = 0;
    int my_rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // keep track of what partitions are my neighbours
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    *send_count = calloc(size, sizeof(int));
    *recv_count = calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        (*send_count)[i] = 0;
        (*recv_count)[i] = 0;
    }

    // count the cells to send
    int e = -1;
    for (int i = 0; i < el_int_loc; ++i) {
        e = local_global_index[i];
        /**
         * see what cells need to be received and sent
         * if a cell has a neighbour, the neighbour has to be received
         * if a cell has a neighbour, the cell has to be sent
         * exclude the external cells
         */
        for (int j = 0; j < 6; ++j) {
            int l = lcc[e][j];
            if (l < el_int_tot) {
                int p = (int) global_local_index[l];
                if (p != my_rank) {
                    ++((*send_count)[p]);
                    ++((*recv_count)[p]);
                }
            }
        }
    }

    // count neighbours
    (*neighbors_count) = 0;
    for (int i = 0; i < size; ++i) {
        if ((*send_count)[i] > 0) {
            ++(*neighbors_count);
        }
    }

    // allocate the lists if no data should be sent set pointer to NULL
    *send_list = (int**) malloc(size * sizeof(int*));
    *recv_list = (int**) malloc(size * sizeof(int*));
    for (int i = 0; i < size; ++i) {
        if ((*send_count)[i] > 0) {
            (*send_list)[i] = (int*) malloc((*send_count)[i] * sizeof(int));
        } else {
            (*send_list)[i] = NULL;
        }
        if ((*recv_count)[i] > 0) {
            (*recv_list)[i] = (int*) malloc((*recv_count)[i] * sizeof(int));
        } else {
            (*recv_list)[i] = NULL;
        }
    }

    // list progr keeps track of the progress of each list
    int* s_list_progr = (int*) calloc(size, sizeof(int));
    int* r_list_progr = (int*) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        s_list_progr[i] = 0;
        r_list_progr[i] = 0;
    }

    // initalize the lists
    for (int i = 0; i < el_int_loc; ++i) {
        e = local_global_index[i];
        for (int j = 0; j < 6; ++j) {
            int l = lcc[e][j];
            if (l < el_int_tot) {
                int p = (int) global_local_index[l];
                if (p != my_rank) {
                    (*send_list)[p][s_list_progr[p]] = i;
                    ++(s_list_progr[p]);
                    (*recv_list)[p][r_list_progr[p]] = l;
                    ++(r_list_progr[p]);
                }
            }
        }
    }

    free(s_list_progr);
    free(r_list_progr);

    return result;
}

/**
 * The external cells are not needed to be taken into accout because their
 * value will always be 0, therefore this can be done directly in the 
 * computation part
 */
int distr_shrink(int* local_global, int el_int_loc, double** array) {
    int result = 0;

    double* tmp = (double*) calloc(el_int_loc, sizeof(double));
    for (int i = 0; i < el_int_loc; ++i) {
        tmp[i] = (*array)[local_global[i]];
    }
    free(*array);
    *array = tmp;

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
                   idx_t** epart, idx_t** npart, idx_t* objval) {
    /********** START INITIALIZATION **********/
    int i = 0;
    int el_int_loc = 0;
    int el_int_tot = 0;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /** Partition data start */
    int part_res = -1;
    // perform partition only on proc 0
    if (my_rank == 0) {
        part_res =  partition(file_in, part_type, nintci, nintcf, nextci, nextcf,
                              lcc, bs, be, bn, bw, bl, bh, bp, su,
                              points_count, points, &el_int_tot, elems,
                              global_local_index, epart, npart, objval);
    }

    MPI_Bcast(&el_int_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(points_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // init arrays before receiving them
    if (my_rank != 0) {
        *global_local_index = (int*) calloc(el_int_tot, sizeof(int));
        *elems = (int*) malloc(el_int_tot * 8 * sizeof(int));
        *lcc  = (int**) malloc(el_int_tot * sizeof(int*));
        **lcc = (int*)  malloc(el_int_tot * 6 * sizeof(int));
        *points  = (int**) malloc(*points_count * sizeof(int*));
        **points = (int*)  malloc(*points_count * 3 * sizeof(int));
    }

    // broadcast the necessary arrays that do not need to be partitioned
    MPI_Bcast(*global_local_index, el_int_tot, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(*elems, el_int_tot * 8, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(**lcc, el_int_tot * 6, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < el_int_tot; ++i) {
        (*lcc)[i] = &((**lcc)[i * 6]);
    }
    MPI_Bcast(**points, *points_count * 3, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < *points_count; ++i) {
        (*points)[i] = &((**points)[i * 3]);
    }
    map_local_global(el_int_tot, *global_local_index,
                     local_global_index, &el_int_loc);

    // initialize the arrays
    *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));
    for (i = 0; i <= 10; i++) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    *var = (double*) calloc((*nextcf) + 1, sizeof(double));
    *cgup = (double*) calloc((*nextcf) + 1, sizeof(double));
    if (my_rank == 0) {
        for (i = 0; i < el_int_tot; ++i) {
            (*cgup)[i] = 0.0;
            (*var)[i] = 0.0;
        }
        for (i = (*nextci); i <= (*nextcf); i++) {
            (*var)[i] = 0.0;
            (*cgup)[i] = 0.0;
            (*bs)[i] = 0.0;
            (*be)[i] = 0.0;
            (*bn)[i] = 0.0;
            (*bw)[i] = 0.0;
            (*bl)[i] = 0.0;
            (*bh)[i] = 0.0;
        }
        for (i = 0; i <= (*nintcf); i++)
            (*cgup)[i] = 1.0 / ((*bp)[i]);
    } else {
        *bs = (double*) calloc((*nextcf) + 1, sizeof(double));
        *be = (double*) calloc((*nextcf) + 1, sizeof(double));
        *bn = (double*) calloc((*nextcf) + 1, sizeof(double));
        *bw = (double*) calloc((*nextcf) + 1, sizeof(double));
        *bl = (double*) calloc((*nextcf) + 1, sizeof(double));
        *bh = (double*) calloc((*nextcf) + 1, sizeof(double));
        *bp = (double*) calloc((*nextcf) + 1, sizeof(double));
        *su = (double*) calloc((*nextcf) + 1, sizeof(double));
    }
    MPI_Bcast(*bs, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*be, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*bn, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*bw, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*bl, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*bh, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*bp, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*su, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*cgup, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*var, (*nextcf) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    distr_shrink(*local_global_index, el_int_loc, bs);
    distr_shrink(*local_global_index, el_int_loc, be);
    distr_shrink(*local_global_index, el_int_loc, bn);
    distr_shrink(*local_global_index, el_int_loc, bw);
    distr_shrink(*local_global_index, el_int_loc, bl);
    distr_shrink(*local_global_index, el_int_loc, bh);
    distr_shrink(*local_global_index, el_int_loc, bp);
    distr_shrink(*local_global_index, el_int_loc, su);
    distr_shrink(*local_global_index, el_int_loc, cgup);
    distr_shrink(*local_global_index, el_int_loc, var);

    /** Partition data end */

    /** set up communication */
    init_commlist(el_int_loc, *local_global_index,
                  el_int_tot, *global_local_index,
                  *lcc,
                  neighbors_count,
                  send_count, send_list,
                  recv_count, recv_list);

    *nintcf = el_int_loc;

    return 0;
}

