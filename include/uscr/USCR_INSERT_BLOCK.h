#ifndef USCR_INSERT_BLOCK_H
#define USCR_INSERT_BLOCK_H
#include "blas_enum.h"
#include "INSERTING.h"
//int BLAS_iuscr_insert_block(blas_sparse_matrix A,SCALAR_IN val ,int i,int j);
//int BLAS_suscr_insert_block(blas_sparse_matrix A,SCALAR_IN val ,int i,int j);
int BLAS_duscr_insert_block(d_matrix* A,double** val ,int row_stride,int clo_stride, int bi,int bj);
//int BLAS_cuscr_insert_block(blas_sparse_matrix A,SCALAR_IN val ,int i,int j);
//int BLAS_zuscr_insert_block(blas_sparse_matrix A,SCALAR_IN val ,int i,int j);
#endif // USCR_INSERT_BLOCK_H
