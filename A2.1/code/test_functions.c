/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include "metis.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      idx_t ne, double *cgup, idx_t* epart, idx_t* npart,
                      idx_t* objval) {
    int metis_result = METIS_ERROR;

    idx_t nn = 0;
    idx_t* eptr = NULL;
    idx_t* eind = NULL;
    idx_t* vwgt = NULL; // can stay NULL
    idx_t* vsize = NULL; // can stay NULL
    idx_t* ncommon = NULL;
    idx_t* nparts = NULL;
    real_t* tpwgts = NULL; // can stay NULL
    idx_t* options = NULL; // can stay NULL

    metis_result = METIS_PartMeshDual(&ne, &nn, eptr, eind, vwgt, vsize, 
                                      ncommon, nparts, tpwgts, options, 
                                      objval, epart, npart);

    return metis_result;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int *num_elems,
                       int neighbors_count, int* send_count, int** send_list, int* recv_count,
                       int** recv_list) {
    // Return an error if not implemented
    return -1;
}

