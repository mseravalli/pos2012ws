/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include "test_functions.h"

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      idx_t ne, double *cgup, idx_t* epart, idx_t* npart,
                      idx_t* objval) {
    // Return an error if not implemented
    return -1;
}

int test_communication(char *file_in, char *file_vtk_out, int *local_global_index, int *num_elems,
                       int neighbors_count, int* send_count, int** send_list, int* recv_count,
                       int** recv_list) {
    // Return an error if not implemented
    return -1;
}

