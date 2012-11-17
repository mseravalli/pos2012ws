/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>

#include "util_read_files.h"
#include "initialization.h"

int partition(char* part_type, 
              int elems_count, int points_count, 
              int* elems,
              idx_t ncommon, idx_t nparts,
              idx_t** objval, idx_t** epart, idx_t** npart) {
    int result = METIS_OK;

    //idx_t ne = elems_count;
    //idx_t nn = points_count;

    //idx_t* eptr = malloc((ne + 2) * sizeof(idx_t));
    //idx_t* eind = malloc((ne + 1) * 8 * sizeof(idx_t));
    //idx_t* vwgt = NULL;
    //idx_t* vsize = NULL;
    //real_t* tpwgts = NULL;
    //idx_t options[METIS_NOPTIONS];

    //*epart = malloc(ne * sizeof(int));
    //*npart = malloc(nn * sizeof(int));

    //// init eprt and eind
    //for (int i = 0; i < ne + 2; ++i) {
    //    eptr[i] = (idx_t) (i * 8);
    //}
    //for (int i = 0; i < (ne + 1) * 8; ++i) {
    //    eptr[i] = (idx_t) elems[i];
    //}

    printf("call the metis function\n");
    
    //METIS_SetDefaultOptions(options);
    
    //result = METIS_PartMeshDual(&ne, &nn,
    //                            eptr, eind,
    //                            vwgt, vsize,
    //                            &ncommon, &nparts,
    //                            tpwgts, options,
    //                            *objval, *epart, *npart);

    switch(result) {
        case METIS_OK:
            printf("partitined correctly\n");
            result = 0;
            break;
        case METIS_ERROR_INPUT:
            printf("Metis input error \n");
            break;
        case METIS_ERROR_MEMORY:
            printf("Metis memory error \n");
            break;
        case METIS_ERROR:
            printf("Metis error \n");
            break;
        default:
            break;
    }
    
    fflush(stdout);

    return result;
}

int initialization(char* file_in, char* part_type, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int* neighbors_count, int** send_count, int*** send_list, int** recv_count,
                   int*** recv_list, idx_t** epart, idx_t** npart, idx_t** objval) {
    /********** START INITIALIZATION **********/
    int i = 0;
    // read-in the input file
    int f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                                   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
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

    idx_t ncommon = 4;
    idx_t nparts = 2;
    int part_result = partition(part_type, 
                               (*nintcf - *nintci), *points_count, 
                               *elems,
                               ncommon, nparts,
                               objval, epart, npart);

    return part_result;
}

