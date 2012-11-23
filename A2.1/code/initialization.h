/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_
#include "metis.h"

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, idx_t** epart, idx_t** npart, idx_t* objval,
                   int* local_elems);

int init_commlist(int local_elems, int* local_global_index, // in, in
                  int elems_count, idx_t* epart, int** lcc, // in, in, in
                  int** commlist);                        // out

#endif /* INITIALIZATION_H_ */

