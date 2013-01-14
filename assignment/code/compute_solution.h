/*
 * compute_solution.h
 *
 *  Created on: Oct 21, 2012
 *      Author: petkovve
 */

#ifndef COMPUTE_SOLUTION_H_
#define COMPUTE_SOLUTION_H_

#define TAG_BP 10
#define TAG_BS 20
#define TAG_BW 30
#define TAG_BL 40
#define TAG_BN 50
#define TAG_BE 60
#define TAG_BH 70

#define TAG_VAR  100
#define TAG_SU   110
#define TAG_CGUP 120

#define TAG_DIR1 200
#define TAG_DIR2 200

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list, int* recv_count, int** recv_list,
                     double* original_b, int** original_lcc);

#endif /* COMPUTE_SOLUTION_H_ */

