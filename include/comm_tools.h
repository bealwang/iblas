#ifndef COMM_TOOLS_H
#define COMM_TOOLS_H
#include "properties.h"
#include <omp.h>
#include <stdio.h>

void* aligned_malloc(size_t size);
void aligned_free(void* p);
void** malloc2d(int rows, int cols, int size);
void initarray(int **array,int rows,int cols);
//p is the input matrix,a1 a2 is the new order,and m n is the size of matrix
void matrix_resort(int** p,int* a1,int* a2,int m,int n);
int minval(int* a,int n_a);
int maxval(int* a,int n_a);
int count(double* a,int n_a,int e);
void dump_matrix(void*);
#endif // COMM_TOOLS_H
