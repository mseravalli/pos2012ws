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

#include "compute_solution.h"

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf,
                     int** lcc, double* bp, double* bs, double* bw, double* bl,
                     double* bn, double* be, double* bh, double* cnorm, 
                     double* var, double *su, double* cgup, 
                     double* residual_ratio, int* local_global_index, 
                     int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list, 
                     int* recv_count, int** recv_list) {
    //TODO: be sure that the following statements do not crash everything
    --nintcf;

    // communication initialization start
    // define the customized types
    int size, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
    struct requests* req = (struct requests*) calloc(size, sizeof(struct requests));
    // communication initialization end

    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;

    // allocate arrays used in gccg
    int nomax = 3;

    /** the reference residual*/
    double resref = 0.0;
    double global_resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    // exchange su
    // communication start
    for (int i = 0; i < size; ++i) {
        if (send_count[i] > 0) {
            MPI_Send(su,   1, send_types[i], i, TAG_SU,   MPI_COMM_WORLD);
        }
        if (recv_count[i] > 0) {
            MPI_Recv(su,   1, recv_types[i], i, TAG_SU,   MPI_COMM_WORLD, &status);
        }
    }
        // communication end
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
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/

        // communication start
        for (int i = 0; i < size; ++i) {
            if (send_count[i] > 0) {
                MPI_Isend(cgup, 1, send_types[i], i, TAG_CGUP, MPI_COMM_WORLD, &(req[i].cgup_send));
                MPI_Isend(bp, 1, send_types[i], i, TAG_BP, MPI_COMM_WORLD, &(req[i].bp_send));
                MPI_Isend(bs, 1, send_types[i], i, TAG_BS, MPI_COMM_WORLD, &(req[i].bs_send));
                MPI_Isend(bw, 1, send_types[i], i, TAG_BW, MPI_COMM_WORLD, &(req[i].bw_send));
                MPI_Isend(bl, 1, send_types[i], i, TAG_BL, MPI_COMM_WORLD, &(req[i].bl_send));
                MPI_Isend(bn, 1, send_types[i], i, TAG_BN, MPI_COMM_WORLD, &(req[i].bn_send));
                MPI_Isend(be, 1, send_types[i], i, TAG_BE, MPI_COMM_WORLD, &(req[i].be_send));
                MPI_Isend(bh, 1, send_types[i], i, TAG_BH, MPI_COMM_WORLD, &(req[i].bh_send));

                MPI_Isend(var,  1, send_types[i], i, TAG_VAR,  MPI_COMM_WORLD, &(req[i].var_send));
                MPI_Isend(su,   1, send_types[i], i, TAG_SU,   MPI_COMM_WORLD, &(req[i].su_send));
            }
            if (recv_count[i] > 0) {
                MPI_Irecv(cgup, 1, recv_types[i], i, TAG_CGUP, MPI_COMM_WORLD, &(req[i].cgup_recv));
                MPI_Irecv(bp, 1, recv_types[i], i, TAG_BP, MPI_COMM_WORLD, &(req[i].bp_recv));
                MPI_Irecv(bs, 1, recv_types[i], i, TAG_BS, MPI_COMM_WORLD, &(req[i].bs_recv));
                MPI_Irecv(bw, 1, recv_types[i], i, TAG_BW, MPI_COMM_WORLD, &(req[i].bw_recv));
                MPI_Irecv(bl, 1, recv_types[i], i, TAG_BL, MPI_COMM_WORLD, &(req[i].bl_recv));
                MPI_Irecv(bn, 1, recv_types[i], i, TAG_BN, MPI_COMM_WORLD, &(req[i].bn_recv));
                MPI_Irecv(be, 1, recv_types[i], i, TAG_BE, MPI_COMM_WORLD, &(req[i].be_recv));
                MPI_Irecv(bh, 1, recv_types[i], i, TAG_BH, MPI_COMM_WORLD, &(req[i].bh_recv));

                MPI_Irecv(var,  1, recv_types[i], i, TAG_VAR,  MPI_COMM_WORLD, &(req[i].var_recv));
//              MPI_Irecv(su,   1, recv_types[i], i, TAG_SU,   MPI_COMM_WORLD, &(req[i].su_recv));
            }
        }
        // communication end

        // update the old values of direc wait for cgup to be received
        for (int i = 0; i < size; ++i) {
            if (recv_count[i] > 0) {
                MPI_Wait(&(req[i].cgup_recv), &status);
            }
        }        
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess (approximation) for direc wait for the bs to be received
        for (int i = 0; i < size; ++i) {
            if (recv_count[i] > 0) {
                MPI_Wait(&(req[i].bp_recv), &status);
                MPI_Wait(&(req[i].bs_recv), &status);
                MPI_Wait(&(req[i].bw_recv), &status);
                MPI_Wait(&(req[i].bl_recv), &status);
                MPI_Wait(&(req[i].bn_recv), &status);
                MPI_Wait(&(req[i].be_recv), &status);
                MPI_Wait(&(req[i].bh_recv), &status);
            }
        }        
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] -
                         bs[nc] * direc1[lcc[nc][0]] -
                         bw[nc] * direc1[lcc[nc][3]] -
                         bl[nc] * direc1[lcc[nc][4]] -
                         bn[nc] * direc1[lcc[nc][2]] -
                         be[nc] * direc1[lcc[nc][1]] -
                         bh[nc] * direc1[lcc[nc][5]];
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

            oc1 = global_occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
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
                oc1 = global_occ / cnorm[1];

                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + adxor2[nc] * direc2[nc];
                }
                // communicate the global occ
                MPI_Allreduce(&occ, &global_occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                oc2 = global_occ / cnorm[2];

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
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
        //wait to receive var
        for (int i = 0; i < size; ++i) {
            if (recv_count[i] > 0) {
                MPI_Wait(&(req[i].var_recv), &status);
            }
        }        
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            var[nc] = var[nc] + omega * direc1[nc];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
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
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

    free(resvec);
    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);

    return iter;
}


