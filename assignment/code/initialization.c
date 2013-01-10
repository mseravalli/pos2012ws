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
              int el_int_glob, int points_count,
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    idx_t ne = el_int_glob;
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
              int el_int_glob, int points_count,
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    idx_t ne = el_int_glob;
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
                        int el_int_glob, int points_count,
                        int* elems,
                        idx_t ncommon, idx_t nparts,
                        idx_t* objval, idx_t** epart, idx_t** npart) {
    int result = 0;

    *epart = malloc(el_int_glob * sizeof(idx_t));
    *npart = malloc(points_count * sizeof(idx_t));

    int elems_per_part = (int) el_int_glob / nparts;

    for (int i = 0; i < el_int_glob; ++i) {
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
              int* el_int_glob, int** elems,
              int** part_elems,
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

    *el_int_glob = (*nintcf - *nintci) + 1;

    if (strcmp(part_type, "dual") == 0 && size > 1) {
        part_result = dual_partition(part_type,
                                     *el_int_glob, *points_count,
                                     *elems,
                                     ncommon, nparts,
                                     objval, epart, npart);
    } else if (strcmp(part_type, "nodal") == 0 && size > 1) {
        part_result = nodal_partition(part_type,
                                      *el_int_glob, *points_count,
                                      *elems,
                                      ncommon, nparts,
                                      objval, epart, npart);
    } else {
        part_result = classical_partition(part_type,
                                          *el_int_glob, *points_count,
                                          *elems,
                                          ncommon, nparts,
                                          objval, epart, npart);
    }

    if (part_result != 0) {
        printf("partition failed\n");
        return -1;
    }

    *part_elems = calloc(*el_int_glob, sizeof(int));
    for (int i = 0; i < *el_int_glob; ++i) {
        (*part_elems)[i] = (*epart)[i];
    }

    free(*epart);
    free(*npart);

    return part_result;
}

int build_global_local(int el_int_glob, int* part_elems,
                       int* recv_count, int** recv_list,
                       int* el_int_loc,
                       int* el_ext_loc,
                       int** global_local) {
    *global_local = (int*) calloc(el_int_glob, sizeof(int));
    
    int my_rank = -1;
    int size = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int count = 0;
    for (int i = 0; i < el_int_glob; ++i) {
        if (part_elems[i] == my_rank) {
            (*global_local)[i] = count;
            ++count;
        } else {
            (*global_local)[i] = -1;
        }
    }
    
    *el_int_loc = count;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < recv_count[i]; ++j) {
            if ((*global_local)[recv_list[i][j]] == -1) {
                (*global_local)[recv_list[i][j]] = count;
                ++count;
            }
        }
    }

    *el_ext_loc = count - *el_int_loc; 

    return 0;
}

int build_local_global(int el_int_glob, int* global_local, int** local_global) {
    int count = 0;
    for (int i = 0; i < el_int_glob; ++i) {
        if (global_local[i] != -1) {
            ++count;
        }
    }

    *local_global = (int*) calloc(count, sizeof(int));
    for (int i = 0; i < el_int_glob; ++i) {
        if (global_local[i] != -1) {
            (*local_global)[global_local[i]] = i;
        }
    }
    return 0;
}

int compare(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

int remove_duplicates(int** array, int num, int* new_num) {
    qsort(*array, num, sizeof(int), compare);

    // look how may duplicates there are
    int dupl = 0;
    int last = (*array)[0];
    for (int i = 1; i < num; ++i) {
        if ((*array)[i] == last) {
            ++dupl;
        } else {
            last = (*array)[i];
        }
    }

    if (dupl == 0)
        return 0;

    int* tmp = (int*) calloc(num - dupl, sizeof(int));
    *new_num = num - dupl;

    dupl = 0;
    last = (*array)[0];
    tmp[dupl] = last;
    for (int i = 1; i < num; ++i) {
        if ((*array)[i] != last) {
            ++dupl;
            last = (*array)[i];
            tmp[dupl] = last;
        } 
    }

    free(*array);
    *array = tmp;

    return 0;
}

/**
 * Initialize the communication lists
 */
int init_commlist(int el_int_glob, int* part_elems, int** lcc,  // i,i,i
                  int* neighbors_count,                        // o
                  int** send_count, int*** send_list,          // o, o
                  int** recv_count, int*** recv_list) {        // o, o
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
    for (int i = 0; i < el_int_glob; ++i) {
        // if the element is not owned my the current processor skip it
        if (part_elems[i] != my_rank) {
            continue;
        }
        /**
         * see what cells need to be received and sent
         * if a cell has a neighbour, the neighbour has to be received
         * if a cell has a neighbour, the cell has to be sent
         * exclude the external cells
         */
        for (int j = 0; j < 6; ++j) {
            int l = lcc[i][j];
            if (l < el_int_glob) {
                int p = (int) part_elems[l];
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
    for (int i = 0; i < el_int_glob; ++i) {
        if (part_elems[i] != my_rank) {
            continue;
        }
        for (int j = 0; j < 6; ++j) {
            int l = lcc[i][j];
            if (l < el_int_glob) {
                int p = (int) part_elems[l];
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

//  // remove duplicates
//      for (int i = 0; i < size; ++i) {
//      if ((*send_count)[i] > 0) {
//          remove_duplicates(&((*send_list)[i]), 
//                            (*send_count)[i], 
//                            &((*send_count)[i]));
//      }
//      if ((*recv_count)[i] > 0) {
//          remove_duplicates(&((*recv_list)[i]), 
//                            (*recv_count)[i],
//                            &((*recv_count)[i]));
//      }
//  }

    return result;
}

/**
 * The last element is 0, this correspond to the external cells, lcc will take
 * this into account 
 */
int distr_shrink(int* local_global, 
                 int el_int_loc, 
                 int el_ext_loc, 
                 double** array) {
    int result = 0;

    double* tmp = (double*) calloc(el_int_loc + el_ext_loc + 1, sizeof(double));
    for (int i = 0; i < el_int_loc + el_ext_loc; ++i) {
        tmp[i] = (*array)[local_global[i]];
    }
    tmp[el_int_loc + el_ext_loc] = 0;
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
    int el_int_glob = 0;
    int el_ext_loc = 0;

    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /** Partition data start */
    int part_res = -1;
    int* part_elems = NULL;
    // perform partition only on proc 0
    if (my_rank == 0) {
        part_res =  partition(file_in, part_type, nintci, nintcf, nextci, nextcf,
                              lcc, bs, be, bn, bw, bl, bh, bp, su,
                              points_count, points, &el_int_glob, elems,
                              &part_elems, epart, npart, objval);
    }

    MPI_Bcast(&el_int_glob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(points_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // init arrays before receiving them
    if (my_rank != 0) {
        part_elems = (int*) calloc(el_int_glob, sizeof(int));
        *elems = (int*) malloc(el_int_glob * 8 * sizeof(int));
        *lcc  = (int**) malloc(el_int_glob * sizeof(int*));
        **lcc = (int*)  malloc(el_int_glob * 6 * sizeof(int));
        *points  = (int**) malloc(*points_count * sizeof(int*));
        **points = (int*)  malloc(*points_count * 3 * sizeof(int));
    }

    // broadcast the necessary arrays that do not need to be partitioned
    MPI_Bcast(part_elems, el_int_glob, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(*elems, el_int_glob * 8, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(**lcc, el_int_glob * 6, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < el_int_glob; ++i) {
        (*lcc)[i] = &((**lcc)[i * 6]);
    }
    MPI_Bcast(**points, *points_count * 3, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < *points_count; ++i) {
        (*points)[i] = &((**points)[i * 3]);
    }

    // set up communication
    init_commlist(el_int_glob, part_elems,
                  *lcc,
                  neighbors_count,
                  send_count, send_list,
                  recv_count, recv_list);
    
    build_global_local(el_int_glob, part_elems,
                       *recv_count, *recv_list,
                       &el_int_loc,
                       &el_ext_loc,
                       global_local_index);
  
    build_local_global(el_int_glob, *global_local_index, local_global_index);
    
    // local values for lcc
    for (int i = 0; i < el_int_loc; ++i) {
        int li = (*local_global_index)[i];
        for (int j = 0; j < 6; ++j) {
            if ((*lcc)[li][j] < el_int_glob) {
                (*lcc)[i][j] = (*global_local_index)[(*lcc)[li][j]];
            } else {
                (*lcc)[i][j] = el_int_loc + el_ext_loc;
            }
        }
    }

    // local values for send and recv list
    for (int i = 0; i < size; ++i) { 
        for (int j = 0; j < (*send_count)[i]; ++j) {
            (*send_list)[i][j] = (*global_local_index)[(*send_list)[i][j]];
        }
        for (int j = 0; j < (*recv_count)[i]; ++j) {
            (*recv_list)[i][j] = (*global_local_index)[(*recv_list)[i][j]];
        }
    }

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
        for (i = 0; i < el_int_glob; ++i) {
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

    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bs);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, be);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bn);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bw);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bl);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bh);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, bp);

    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, su);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, cgup);
    distr_shrink(*local_global_index, el_int_loc, el_ext_loc, var);

    /** Partition data end */

    *nintci = 0;
    *nintcf = el_int_loc;
    *nextcf = el_int_loc + el_ext_loc;

    free(part_elems);

    return 0;
}

