/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_
#include "metis.h"

#define RECV_ELEM  5
#define SEND_ELEM  10
#define INNER_ELEM 15

// MPI TAGS


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
                   idx_t** epart, idx_t** npart, idx_t* objval, int* local_elems);

int init_commlist(int local_elems, int* local_global_index, // in, in
                  int elems_count, int* global_local_index, // in, in
                  int** lcc,                                // in
                  int** commlist, int* neighbors_count,     // out, out
                  int** send_count, int** recv_count);      // out, out

#endif /* INITIALIZATION_H_ */

