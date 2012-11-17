/**
 * Initialization step - parse the input file, compute data distribution, 
 * initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdlib.h>
#include "metis.h"

#include "util_read_files.h"
#include "initialization.h"

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf,
                   int* nextci, int* nextcf, int*** lcc, double** bs,
                   double** be, double** bn, double** bw, double** bl,
                   double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup,
                   double** oc, double** cnorm, int** local_global_index,
                   int** global_local_index, int* neighbors_count,
                   int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, idx_t** epart, idx_t** npart, idx_t** objval) {
    /********** START INITIALIZATION **********/
    int i = 0;
    // read-in the input file
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci,
                                   &*nextcf, &*lcc, &*bs, &*be, &*bn, &*bw, &*bl,
                                   &*bh, &*bp, &*su, &*points_count,
                                   &*points, &*elems);

    if ( f_status != 0 ) return f_status;

    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *oc = (double*) calloc(sizeof(double), (*nintcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*oc)[i] = 0.0;
        (*cnorm)[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*cgup)[i] = 0.0;
        (*var)[i] = 0.0;
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bl)[i] = 0.0;
        (*bh)[i] = 0.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ )
        (*cgup)[i] = 1.0 / ((*bp)[i]);

    int metis_result = METIS_ERROR;

    
    idx_t ne = (idx_t)(*nextcf);
    idx_t nn = (idx_t)(*points_count);
    idx_t* eptr = calloc(ne + 1, sizeof(idx_t));
    idx_t* eind = calloc(nn, sizeof(idx_t));
    idx_t* vwgt = NULL; // can stay NULL
    idx_t* vsize = NULL; // can stay NULL
    idx_t ncommon = 1;
    idx_t nparts = 6;
    real_t* tpwgts = NULL; // can stay NULL
    idx_t* options = NULL; // can stay NULL

    metis_result = METIS_PartMeshDual(&ne, &nn, 
                                      eptr, eind, 
                                      vwgt, vsize, 
                                      &ncommon, &nparts, 
                                      tpwgts, options, 
                                      *objval, *epart, *npart);


    return 0;
}

