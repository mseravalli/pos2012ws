#ifndef XWRITE_H_
#define XWRITE_H_

int write_result( char *inFileName, 
                  char *outFileName, 
                  int NINTCI, 
                  int NINTCF, 
                  double *VAR, 
                  int ITER, 
                  double RATIO );

int write_result_vtk( char* outFileName, 
                      int startIdx, 
                      int endIdx, 
                      int nodeCnt, 
                      int** points, 
                      int** elems, 
                      double* vector ); 

int write_perf_data( double* perf_data );

#endif /* XWRITE_H_ */
