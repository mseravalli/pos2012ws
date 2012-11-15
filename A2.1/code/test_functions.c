/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include "parmetis.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      long ne, double *cgup, long* epart, long* npart,
                      long* objval) {
    int metis_result = METIS_ERROR;

    long nn = 0;
    long* eptr = NULL;
    long* eind = NULL;
    long* vwgt = NULL; // can stay NULL
    long* vsize = NULL; // can stay NULL
    long* ncommon = NULL;
    long* nparts = NULL;
    float* tpwgts = NULL; // can stay NULL
    long* options = NULL; // can stay NULL

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

