#ifndef DENSE_H
#define DENSE_H
#include "comm_tools.h"
void dblock_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr);
void dblock_Z_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr);
void dblock_T_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr);
void dblock_H_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr);
void dblock_r_mult_vec (double *A, double *x, int n, double *y, int m, int *ierr);
void dblock_l_mult_vec (double *A, double *x, int n, double *y, int m, int *ierr);
#endif // DENSE_H


