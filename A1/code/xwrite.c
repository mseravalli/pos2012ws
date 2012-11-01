#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xwrite.h"

/**
 * Export the results to a file
 *
 * @param inFileName
 * @param outFileName
 * @param NINTCI
 * @param NINTCF
 * @param VAR
 * @param ITER
 * @param RATIO
 * @return
 */
int write_result( char *inFileName, 
                  char *outFileName, 
                  int NINTCI, 
                  int NINTCF,
                  double *VAR, 
                  int ITER, 
                  double RATIO ) {
    double *IPOINT = ( double * ) malloc( ( NINTCF + 1 ) * sizeof( double ) );
    int I1, I2, I3, I4, I5;

    int NC;
    for ( NC = NINTCI; NC <= NINTCF; NC++ ) {
        IPOINT[NC] = NC;
    }

    if ( NINTCF <= 1 ) {
        printf( "Error: NINTCF <= 1\n" );
        return -1;
    }

    I1 = NINTCF + 1;

    while ( I1 != 0 ) {
        I1 = I1 / 2;
        I2 = NINTCF + 1 - I1;
        I4 = 1;

        do {
            I3 = I4;
            do {
                I5 = I3 + I1;
                if ( VAR[I3] <= VAR[I5] ) {
                    break;
                }

                double ZDUM = VAR[I3], IDUM  = IPOINT[I3];

                VAR[I3] = VAR[I5];
                IPOINT[I3] = IPOINT[I5];
                VAR[I5] = ZDUM;
                IPOINT[I5] = IDUM;
                I3 = I3 - I1;
            } while ( I3 >= 1 );
            I4++;
        } while ( I4 < I2 );
    }

    char out_file_name[128];
    strcpy( out_file_name, outFileName); 
    strcat( out_file_name, ".dat" );
    FILE *fp = fopen( out_file_name, "w" );
    if ( fp == NULL ) {
        printf( "Error opening file %s for writing\n", outFileName );
        return -1;
    }

    printf( "Writing data to %s: ", outFileName );
    fprintf( fp, "               ----------------------------------------\n" );
    fprintf( fp, "               - AVL -  Linear Equation Solver - GCCG -\n" );
    fprintf( fp, "               ----------------------------------------\n" );
    fprintf( fp, "                     diagonal scaling\n" );
    fprintf( fp, "\n\n" );
    fprintf( fp, "     Input File:  %s\n", inFileName );
    fprintf( fp, "     ===========\n\n" );
    fprintf( fp, "     Output File:  %s\n", outFileName );
    fprintf( fp, "     ============\n\n" );
    fprintf( fp, "     No. of Active Cells:  %d\n", NINTCF );
    fprintf( fp, "     ====================\n\n" );
    fprintf( fp, "     Iteration Count - Residual Ratio\n" );
    fprintf( fp, "     ================================\n\n" );
    fprintf( fp, "%d %e\n", ITER, RATIO );
    fprintf( fp, "     Addresses Solution ( Minima )           Addresses Solution ( Maxima )\n" );
    fprintf( fp, "     ===========================           ===========================\n\n" );

    int N;
    for ( N = 1; N <= 10; N++ ) {
        fprintf( fp, 
                 "%lf %lf %lf %lf\n", 
                 IPOINT[N], VAR[N],    
                 IPOINT[NINTCF - N + 1], 
                 VAR[NINTCF - N + 1] );
    }

    fclose( fp );
    printf( "done!\n" );

    free( IPOINT );
    return 0;
}

int write_result_vtk( char* current_sim,
                      char *outFileName, 
                      int startIdx, 
                      int endIdx, 
                      int nodeCnt, 
                      int **points, 
                      int **elems, 
                      double *vector ) {
    int i,j;
    int cellCnt = endIdx - startIdx + 1;
    char vtk_path[128];
    strcpy( vtk_path, current_sim );
    strcat( vtk_path, "." );
    strcat( vtk_path, outFileName );
    FILE *fp = fopen( vtk_path, "w" );
    if( fp == NULL ) {
        printf( "Error opening file %s for writing\n", outFileName );
        return -1;
    }

    fprintf( fp, "# vtk DataFile Version 3.0\n" );// file version and identifier
    fprintf( fp,"vtk output\n" );                 // header
    fprintf( fp,"ASCII\n" );                      // file format
    fprintf( fp,"DATASET UNSTRUCTURED_GRID\n" );  // dataset structure

    fprintf( fp,"POINTS %d float\n",nodeCnt );

    printf( "nodeCnt = %d, startInd= %d, endInd = %d, cellCnt=%d\n",
             nodeCnt, 
             startIdx, 
             endIdx, 
             cellCnt );
    for ( i = 1; i<=nodeCnt/2; i++ ) {
        fprintf( fp, 
                 "%d %d %d ", 
                 points[0][2*i],
                 points[1][2*i],
                 points[2][2*i] );
        fprintf( fp,
                 "%d %d %d\n",
                 points[0][2*i+1],
                 points[1][2*i+1],
                 points[2][2*i+1] );
    }

    fprintf( fp, "\nCELLS %d %d\n", cellCnt, 9*cellCnt );
    for ( i = startIdx; i<= endIdx; i++ ) {
        fprintf( fp,
                 "8 %d %d %d %d %d %d %d %d\n",
                 elems[0][i],
                 elems[1][i],
                 elems[2][i],
                 elems[3][i],
                 elems[4][i],
                 elems[5][i],
                 elems[6][i],
                 elems[7][i] );
    }

    fprintf( fp,"\nCELL_TYPES %d\n",cellCnt );
    for ( i = startIdx; i<= endIdx; i++ ) {
        fprintf( fp,"11\n" );//cell type vtk voxel
    }

    fprintf( fp,"\nCELL_DATA %d\n",cellCnt );
    fprintf( fp,"SCALARS scalars float\nLOOKUP_TABLE default\n" );

    for ( i = startIdx; i<= endIdx; i++ ) {
        fprintf( fp,"%f ",vector[i] );
    }

    fclose( fp );
    return 0;
}


int write_perf_data( char* current_sim,
                     char* format,
                     char* phase,
                     long_long exec_time, 
                     long_long* perf_data ) {

    char stats_path[128];
    strcpy( stats_path, current_sim );
//  strcat( stats_path, "." );
//  strcat( stats_path, format );
    strcat( stats_path, ".pstats.dat" );

    // create a new file every run
    FILE* stats_file = NULL; 
    
    if ( strcmp( phase, "INPUT" ) == 0 ) {
        stats_file = fopen( stats_path, "w" );
    } else {
        stats_file = fopen( stats_path, "a" );
    }

    if( stats_file == NULL ) {
        printf( "Error opening file %s for writing\n", stats_path );
        return -1;
    }

    fprintf( stats_file, 
        "%s_Exec_time=%.4lfs\n", 
        phase, 
        (double) exec_time * 1e-6 );

    fprintf( stats_file, "%s_PAPI_L2_TCM=%ld\n", phase, perf_data[0] );
    fprintf( stats_file, "%s_PAPI_L2_TCA=%ld\n", phase, perf_data[1] );
    double miss_rate = ( double ) perf_data[0] / ( double ) perf_data[1]; 
    fprintf( stats_file, "%s_L2_MissRate=%.2lf%%\n", phase,  100 * miss_rate );

    fprintf( stats_file, "%s_PAPI_FPI_INS=%ld\n", phase, perf_data[2] );
    fprintf( stats_file, "%s_PAPI_TOT_CYC=%ld\n", phase, perf_data[3] );

    double m_flops = ( double ) perf_data[2] / ( double ) exec_time;
    fprintf( stats_file, "%s_Mflops=%.4lf\n", phase, m_flops );

    double util = ( double ) perf_data[3] / ( ( double ) exec_time * 20.0 ); 
    fprintf( stats_file, "%s_Util=%.2lf%%\n", phase, util );

    fclose( stats_file );

    return 0;

}


