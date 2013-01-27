/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "scorep/SCOREP_User.h"


#include "./compute_solution.h"

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf,
                     int** lcc, double* bp, double* bs, double* bw, double* bl,
                     double* bn, double* be, double* bh, double* cnorm,
                     double* var, double *su, double* cgup,
                     double* residual_ratio, int* local_global_index,
                     int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list,
                     int* recv_count, int** recv_list) {
    --nintcf;

    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;

    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // allocate arrays used in gccg
    int nomax = 3;

    /** the reference residual*/
    double resref = 0.0;
    double global_resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    MPI_Allreduce(&resref, &global_resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref = global_resref;

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    /** the computation vectors */
    double *direc1_curr = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc1_next = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    for (int i = 0; i <= nextcf; ++i) {
        direc1_curr[i] = 0;
        direc2[i] = 0;
    }

    // communication initialization start
    // define the customized types
    int* block_len = NULL;
    MPI_Datatype* send_types =
        (MPI_Datatype*) calloc(size, sizeof(MPI_Datatype));
    MPI_Datatype* recv_types =
        (MPI_Datatype*) calloc(size, sizeof(MPI_Datatype));
    for (int i = 0; i < size; ++i) {
        if (send_count[i] > 0) {
            block_len = (int*) calloc(send_count[i], sizeof(int));
            for (int j = 0; j < send_count[i]; ++j) {
                block_len[j] = 1;
            }
            MPI_Type_indexed(send_count[i],
                             block_len,
                             send_list[i],
                             MPI_DOUBLE,
                             &(send_types[i]));
            MPI_Type_commit(&(send_types[i]));
            free(block_len);
        }
        if (recv_count[i] > 0) {
            block_len = (int*) calloc(recv_count[i], sizeof(int));
            for (int j = 0; j < recv_count[i]; ++j) {
                block_len[j] = 1;
            }
            MPI_Type_indexed(recv_count[i],
                             block_len,
                             recv_list[i],
                             MPI_DOUBLE,
                             &(recv_types[i]));
            MPI_Type_commit(&(recv_types[i]));
            free(block_len);
        }
    }
    // variables used for synchronization
    MPI_Status status;
    MPI_Request* s_req = (MPI_Request*) calloc(size, sizeof(MPI_Request));
    MPI_Request* r_req = (MPI_Request*) calloc(size, sizeof(MPI_Request));
    // communication initialization end

    SCOREP_USER_REGION_DEFINE(OA_Phase);
    SCOREP_USER_OA_PHASE_BEGIN(OA_Phase, "OA_Phase", SCOREP_USER_REGION_TYPE_COMMON);

    // update the old values of direc
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        direc1_curr[nc] = direc1_curr[nc] + resvec[nc] * cgup[nc];
    }
    for (int i = 0; i < size; ++i) {
        if (send_count[i] > 0) {
            MPI_Isend(direc1_curr, 1, send_types[i], i, TAG_DIR1, MPI_COMM_WORLD, &(s_req[i]));
        }
    }

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/

        // start receiving
        for (int i = 0; i < size; ++i) {
            if (recv_count[i] > 0) {
                MPI_Irecv(direc1_curr, 1, recv_types[i], i, TAG_DIR1, MPI_COMM_WORLD, &(r_req[i]));
            }
        }
        for (int i = 0; i < size; ++i) {
            if (send_count[i] > 0) {
                MPI_Wait(&(s_req[i]), &status);
            }
        }
        for (int i = 0; i < size; ++i) {
            if (recv_count[i] > 0) {
                MPI_Wait(&(r_req[i]), &status);
            }
        }
        // end receiving

        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1_curr[nc] -
                         bs[nc] * direc1_curr[lcc[nc][0]] -
                         bw[nc] * direc1_curr[lcc[nc][3]] -
                         bl[nc] * direc1_curr[lcc[nc][4]] -
                         bn[nc] * direc1_curr[lcc[nc][2]] -
                         be[nc] * direc1_curr[lcc[nc][1]] -
                         bh[nc] * direc1_curr[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        double global_occ = 0;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }

            // communicate the global occ
            MPI_Allreduce(&occ,
                          &global_occ,
                          1,
                          MPI_DOUBLE,
                          MPI_SUM,
                          MPI_COMM_WORLD);
            occ = global_occ;
            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1_curr[nc] = direc1_curr[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + adxor1[nc] * direc2[nc];
                }

                // communicate the global occ
                MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                occ = global_occ;
                oc1 = occ / cnorm[1];

                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + adxor2[nc] * direc2[nc];
                }
                // communicate the global occ
                MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                occ  = global_occ;
                oc2 = occ / cnorm[2];

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                    direc1_curr[nc] = direc1_curr[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double global_norm = 0;
        double omega = 0;
        double global_omega = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        // communicate the global omega and cnorm
        MPI_Allreduce(&omega, &global_omega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        omega = global_omega;
        MPI_Allreduce(&(cnorm[nor]), &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        cnorm[nor] = global_norm;

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        double global_res = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            var[nc] = var[nc] + omega * direc1_curr[nc];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
        }
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1_next[nc] = direc1_curr[nc] + resvec[nc] * cgup[nc];
        }
        for (int i = 0; i < size; ++i) {
            if (send_count[i] > 0) {
                MPI_Isend(direc1_next, 1, send_types[i], i, TAG_DIR1, MPI_COMM_WORLD, &(s_req[i]));
            }
        }

        // communicate the residual
        MPI_Allreduce(&res_updated,
                      &global_res,
                      1,
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);
        res_updated = global_res;

        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1_curr[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1_curr[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
        double* tmp = direc1_next;
        direc1_next = direc1_curr;
        direc1_curr = tmp;
    }

    // start receiving
    for (int i = 0; i < size; ++i) {
        if (recv_count[i] > 0) {
            MPI_Irecv(direc1_curr, 1, recv_types[i], i, TAG_DIR1, MPI_COMM_WORLD, &(r_req[i]));
        }
    }
    for (int i = 0; i < size; ++i) {
        if (send_count[i] > 0) {
            MPI_Wait(&(s_req[i]), &status);
        }
    }
    for (int i = 0; i < size; ++i) {
        if (recv_count[i] > 0) {
            MPI_Wait(&(r_req[i]), &status);
        }
    }
    // end receiving

    SCOREP_USER_OA_PHASE_END(OA_Phase);

    free(resvec);
    free(direc1_curr);
    free(direc1_next);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);

    return iter;
}


