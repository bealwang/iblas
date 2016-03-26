#ifndef MIC_USMM_H
#define MIC_USMM_H
#include "comm_tools.h"
#include "blas_enum.h"
#include "types.h"
#include "mic_mbv.h"
//int BLAS_iusmm(enum blas_order_type order,enum blas_trans_type transa,int nrhs,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY B,int ldb,ARRAY C,int ldc);

//int BLAS_susmm(enum blas_order_type order,enum blas_trans_type transa,int nrhs,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY B,int ldb,ARRAY C,int ldc);

int BLAS_dusmm(blas_order_type order, blas_trans_type transa,int nrhs,double alpha,DSPMAT* A, double *B,int ldb,double *C,int ldc);

//int BLAS_cusmm(enum blas_order_type order,enum blas_trans_type transa,int nrhs,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY B,int ldb,ARRAY C,int ldc);

//int BLAS_zusmm(enum blas_order_type order,enum blas_trans_type transa,int nrhs,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY B,int ldb,ARRAY C,int ldc);
#endif // USMM_H
