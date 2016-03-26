#ifndef USCR_BEGIN_H
#define USCR_BEGIN_H
#include "blas_enum.h"
#include "INSERTING.h"

i_matrix* BLAS_iuscr_begin(int m,int n);
d_matrix* BLAS_duscr_begin(int m,int n);
d_matrix* BLAS_duscr_block_begin(int Mb,int Nb,int k,int l);
d_matrix* BLAS_duscr_variable_block_begin(int Mb,int Nb,const int* k,const int* l);
//blas_sparse_matrix BLAS_cuscr_begin(int m,int n);
//blas_sparse_matrix BLAS_zuscr_begin(int m,int n);
#endif // USCR_BEGIN_H
