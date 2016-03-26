#ifndef USCR_INSERT_CLIQUE_H
#define USCR_INSERT_CLIQUE_H
#include "blas_enum.h"
#include "INSERTING.h"
//int BLAS_iuscr_insert_clique(blas_sparse_matrix A,const int k,const int l,const ARRAY val,const int row_stride,const int col_stride,const int *indx,const int *jndx);

//int BLAS_suscr_insert_clique(blas_sparse_matrix A,const int k,const int l,const ARRAY val,const int row_stride,const int col_stride,const int *indx,const int *jndx);

int BLAS_duscr_insert_clique(d_matrix* A,const int k,const int l,const double** val,const int row_stride,const int col_stride,const int* indx,const int* jndx);

//int BLAS_cuscr_insert_clique(blas_sparse_matrix A,const int k,const int l,const ARRAY val,const int row_stride,const int col_stride,const int *indx,const int *jndx);

//int BLAS_zuscr_insert_clique(blas_sparse_matrix A,const int k,const int l,const ARRAY val,const int row_stride,const int col_stride,const int *indx,const int *jndx);
#endif // USCR_INSERT_CLIQUE_H
