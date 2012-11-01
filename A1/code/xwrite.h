#ifndef XWRITE_H_
#define XWRITE_H_

#include <papi.h>

int write_result( char *inFileName, 
                  char *outFileName, 
                  int NINTCI, 
                  int NINTCF, 
                  double *VAR, 
                  int ITER, 
                  double RATIO );

int write_result_vtk( char* current_sim,
                      char* outFileName, 
                      int startIdx, 
                      int endIdx, 
                      int nodeCnt, 
                      int** points, 
                      int** elems, 
                      double* vector ); 

int write_perf_data( char* current_sim,
                     char* format,
                     char* phase, 
                     long_long exec_time, 
                     long_long* perf_data );

#endif /* XWRITE_H_ */
