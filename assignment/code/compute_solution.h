/*
 * compute_solution.h
 *
 *  Created on: Oct 21, 2012
 *      Author: petkovve
 */

#ifndef COMPUTE_SOLUTION_H_
#define COMPUTE_SOLUTION_H_

#include <mpi.h>

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

struct requests {
    MPI_Request bp_send;
    MPI_Request bs_send;
    MPI_Request bw_send;
    MPI_Request bl_send;
    MPI_Request bn_send;
    MPI_Request be_send;
    MPI_Request bh_send;
    MPI_Request var_send;
    MPI_Request su_send;
    MPI_Request cgup_send;

    MPI_Request bp_recv;
    MPI_Request bs_recv;
    MPI_Request bw_recv;
    MPI_Request bl_recv;
    MPI_Request bn_recv;
    MPI_Request be_recv;
    MPI_Request bh_recv;
    MPI_Request var_recv;
    MPI_Request su_recv;
    MPI_Request cgup_recv;
};

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list, int* recv_count, int** recv_list);

#endif /* COMPUTE_SOLUTION_H_ */

