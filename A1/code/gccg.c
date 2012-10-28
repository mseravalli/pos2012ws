#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <papi.h>

#include "vol2mesh.h"
#include "xread.h"
#include "xwrite.h"

#define NUM_EVENTS  6 

/**
 *
 */
int main( int argc, char *argv[] ) {
    if ( argc < 4 )    {
        printf("Usage: %s format input_file output_file\n", argv[0] );
        return EXIT_FAILURE;
    }

    enum FORMAT_T file_format;

    if ( strcmp( argv[1], "bin" ) == 0 ) {
        file_format = FORMAT_BINARY; 
    } else if ( strcmp( argv[1], "text" ) == 0 ) {
        file_format = FORMAT_TEXT; 
    } else {
        printf("Format '%s' not valid please use 'bin' or 'text'\n", argv[1] );
        return EXIT_FAILURE;
    }

    char *file_in = argv[2];
    char *file_out = argv[3];

    // PAPI initialisation start

    int Events[NUM_EVENTS] = { PAPI_L2_TCM, 
                               PAPI_L2_TCA, 
                               PAPI_FP_INS,
                               PAPI_TOT_CYC };
    long_long values[NUM_EVENTS];

    long_long start_usec, end_usec;

    if ( PAPI_library_init( PAPI_VER_CURRENT ) != PAPI_VER_CURRENT ) {
        printf("Error while initialising PAPI\n" );
        return EXIT_FAILURE;
    }

    // PAPI initialisation end

    int status = 0;

    // internal cells start and end index
    int nintci, nintcf;
    // external cells start and end index. The external cells are only ghost 
    // cells. They are accessed only through internal cells
    int nextci, nextcf;
    // link cell-to-cell array. Stores topology information
    int **lcc;
    // red-black colouring of the cells
    int *nboard;

    // boundary coefficients for each volume cell
    double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;

    // gc initialization 
    if ( PAPI_start_counters( Events, 4 ) != PAPI_OK ) {
        printf("Error while using hardware counters\n" );
        return EXIT_FAILURE;
    }
     
    start_usec = PAPI_get_real_usec();

    // read-in the input file
    int f_status = -1;
    if ( file_format == FORMAT_BINARY ) {
        f_status = read_binary( file_in, 
                                &nintci, 
                                &nintcf, 
                                &nextci, 
                                &nextcf, 
                                &lcc,
                                &bs, 
                                &be, 
                                &bn, 
                                &bw, 
                                &bl, 
                                &bh, 
                                &bp, 
                                &su, 
                                &nboard );
    }
    if ( file_format == FORMAT_TEXT ) {
        f_status = read_formatted( file_in, 
                                   &nintci, 
                                   &nintcf, 
                                   &nextci, 
                                   &nextcf, 
                                   &lcc,
                                   &bs, 
                                   &be, 
                                   &bn, 
                                   &bw, 
                                   &bl, 
                                   &bh, 
                                   &bp, 
                                   &su, 
                                   &nboard );
    }

    if ( f_status != 0 ) {
        printf( "failed to initialize data!\n" );
        return EXIT_FAILURE;
    }

    // allocate arrays used in gccg
    int nomax = 3;
    // the reference residual
    double resref = 0.0;
    // the ratio between the reference and the current residual
    double ratio;

    // array storing residuals
    double* resvec = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    // the variation vector -> keeps the result in the end 
    double* var = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );

    // gc* the computation vectors 
    double* direc1 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double* direc2 = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );

    // gc* additional vectors 
    double* cgup = ( double * ) calloc( sizeof( double ), ( nextcf + 1 ) );
    double* oc = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double* cnorm = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double* adxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double* adxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double* dxor1 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );
    double* dxor2 = ( double * ) calloc( sizeof( double ), ( nintcf + 1 ) );

    // initialize the reference residual
    for ( int nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }
    resref = sqrt( resref );
    if ( resref < 1.0e-15 ) {
        printf( "i/o - error: residue sum less than 1.e-15 - %lf\n", resref );
        return EXIT_FAILURE;
    }

    // initialize the arrays
    for ( int nc = 0; nc <= 10; nc++ ) {
        oc[nc] = 0.0;
        cnorm[nc] = 1.0;
    }

    for ( int nc = nintci; nc <= nintcf; nc++ ) {
        cgup[nc] = 0.0;
        var[nc] = 0.0;
    }

    for ( int nc = nextci; nc <= nextcf; nc++ ) {
        var[nc] = 0.0;
        cgup[nc] = 0.0;
        direc1[nc] = 0.0;
        bs[nc] = 0.0;
        be[nc] = 0.0;
        bn[nc] = 0.0;
        bw[nc] = 0.0;
        bl[nc] = 0.0;
        bh[nc] = 0.0;
    }

    for ( int nc = nintci; nc <= nintcf; nc++ ) {
        cgup[nc] = 1.0 / bp[nc];
    }

    int if1 = 0;
    int if2 = 0;
    int iter = 1;
    int nor = 1;
    int nor1 = nor - 1;
    //gc finished initalization 

    if ( PAPI_read_counters( values, NUM_EVENTS ) != PAPI_OK ) {
        printf("Error while using hardware counters\n" );
        return EXIT_FAILURE;
    }

    end_usec = PAPI_get_real_usec();

    if ( write_perf_data ( file_out, 
                           argv[1],
                           "INPUT", 
                           end_usec - start_usec, 
                           values ) != 0 ) {
        printf( "error when trying to write performance data\n" );
    }

    start_usec = PAPI_get_real_usec();

    //gc start computation loop 
    while ( iter < 10000 ) {

        //gc start phase 1 

        // update the old values of direc
        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess ( approximation ) for direc
        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
                - bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
                - bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
                - bh[nc] * direc1[lcc[5][nc]];
        } // gc end phase 1 

        // gc  start phase 2
        // execute normalization steps
        double oc1, oc2, occ;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }
            oc1 = occ / cnorm[1];
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }
            if1++;

        } else if ( nor1 == 2 ) {
            oc1 = 0;
            occ = 0;
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor1[nc] * direc2[nc];
            }

            oc1 = occ / cnorm[1];
            oc2 = 0;
            occ = 0;
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + adxor2[nc] * direc2[nc];
            }

            oc2 = occ / cnorm[2];
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
            }

            if2++;
        }

        cnorm[nor] = 0;
        double omega = 0;

        // compute the new residual
        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        omega = omega / cnorm[nor];

        double resnew = 0.0;
        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            var[nc] = var[nc] + omega * direc1[nc];
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            resnew = resnew + resvec[nc] * resvec[nc];
        }
        resnew = sqrt( resnew );
        ratio = resnew / resref;

        // exit on no improvements of residual
        if ( ratio <= 1.0e-10 ) {
            break;
        }

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( int nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }

            } else if ( nor == 2 ) {
                for ( int nc = nintci; nc <= nintcf; nc++ ) {
                    dxor2[nc] = direc1[nc];
                    adxor2[nc] = direc2[nc];
                }
            }
            nor++;
        }
        nor1 = nor - 1;

    } // gc end phase 2 

    // gc finished computation loop 

    if ( PAPI_read_counters( values, NUM_EVENTS ) != PAPI_OK ) {
        printf("Error while using hardware counters\n" );
        return EXIT_FAILURE;
    }

    end_usec = PAPI_get_real_usec();

    if ( write_perf_data ( file_out, 
                           argv[1],
                           "CALC", 
                           end_usec - start_usec, 
                           values ) != 0 ) {
        printf( "error when trying to write performance data\n" );
    }

    start_usec = PAPI_get_real_usec();

    // gc write output file  
    if ( write_result( file_in, 
                       file_out, 
                       nintci, 
                       nintcf, 
                       var, 
                       iter, 
                       ratio ) != 0 ) {
        printf( "error when trying to write to file %s\n", file_out );
    }

    // convert to mesh
    int** points;
    int** elems;
    int node_cnt;

    vol2mesh( nintci, 
              nintcf, 
              lcc, 
              &node_cnt, 
              &points, 
              &elems);

    file_out = "SU.vtk";
    if ( write_result_vtk( file_out, 
                           nintci, 
                           nintcf, 
                           node_cnt,
                           points, 
                           elems, 
                           su ) != 0 ) {
        printf( "error when trying to write to file %s\n", file_out );
    }

    file_out = "VAR.vtk";
    if ( write_result_vtk( file_out, 
                           nintci, 
                           nintcf, 
                           node_cnt,
                           points, 
                           elems, 
                           var ) != 0 ) {
        printf( "error when trying to write to file %s\n", file_out );
    }

    file_out = "CGUP.vtk";
    if ( write_result_vtk( file_out, 
                           nintci, 
                           nintcf, 
                           node_cnt,
                           points, 
                           elems, 
                           cgup ) != 0 ) {
        printf( "error when trying to write to file %s\n", file_out );
    }

    if ( PAPI_stop_counters( values, NUM_EVENTS ) != PAPI_OK ) {
        printf("Error while using hardware counters\n" );
        return EXIT_FAILURE;
    }

    end_usec = PAPI_get_real_usec();

    if ( write_perf_data ( argv[3],
                           argv[1],
                           "OUTPUT", 
                           end_usec - start_usec, 
                           values ) != 0 ) {
        printf( "error when trying to write performance data\n" );
    }

    //gc Free all the dynamically allocated memory
    free( direc2 ); 
    free( direc1 ); 
    free( dxor2 ); 
    free( dxor1 ); 
    free( adxor2 ); 
    free( adxor1 );
    free( cnorm ); 
    free( oc ); 
    free( var ); 
    free( cgup ); 
    free( resvec ); 
    free( su ); 
    free( bp );
    free( bh ); 
    free( bl ); 
    free( bw ); 
    free( bn ); 
    free( be ); 
    free( bs );

    printf( "Simulation completed successfully!\n" );
    return EXIT_SUCCESS;
}
