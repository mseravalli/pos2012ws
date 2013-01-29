/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef FINALIZATION_H_
#define FINALIZATION_H_

#define TAG_LEN 200
#define TAG_ARR 210
#define TAG_LG  220

/**
 * The function provides the complete array only to the master 
 */
int merge_arrayd(int length_loc, double* array_loc, int* local_global, 
                 double** array_glob);

/**
 * Runs mainly on the master process
 */
void finalization(char* file_in, char* out_prefix, int total_iters, double residual_ratio,
                  int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
                  double* cgup, double* su, int* local_global);

#endif /* FINALIZATION_H_ */

