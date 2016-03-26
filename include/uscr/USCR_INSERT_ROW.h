#ifndef USCR_INSERT_ROW_H
#define USCR_INSERT_ROW_H
#include "blas_enum.h"
#include "INSERTING.h"
//int BLAS_iuscr_insert_row(blas_sparse_matrix A,int i,int nz,const ARRAY val,const int *indx);
//int BLAS_suscr_insert_row(blas_sparse_matrix A,int i,int nz,const ARRAY val,const int *indx);
int BLAS_duscr_insert_row(d_matrix*,int i,int nz,const double* val,const int* jndx);
//int BLAS_cuscr_insert_row(blas_sparse_matrix A,int i,int nz,const ARRAY val,const int *indx);
//int BLAS_zuscr_insert_row(blas_sparse_matrix A,int i,int nz,const ARRAY val,const int *indx);
#endif // USCR_INSERT_ROW_H
