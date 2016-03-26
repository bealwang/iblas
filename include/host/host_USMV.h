#ifndef HOST_USMV_H
#define HOST_USMV_H
#include "blas_enum.h"
#include "host_mbv.h"
#include "link.h"
//int BLAS_iusmv(enum blas_trans_type transa,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY x,int incx, ARRAY y,int incy);
//int BLAS_susmv(enum blas_trans_type transa,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY x,int incx, ARRAY y,int incy);
int BLAS_dusmv(blas_trans_type transa,double alpha,dsp_linknode* A, double* x, int incx, double* y, int incy);
//int BLAS_dusmv (dsp_linknode* A, const double* x, int incx, double* y, int incy,blas_trans_type transa=ORIGIN_MATRIX,int alpha=1);
//int BLAS_cusmv(enum blas_trans_type transa,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY x,int incx, ARRAY y,int incy);
//int BLAS_zusmv(enum blas_trans_type transa,SCALAR_IN alpha,blas_sparse_matrix A,const ARRAY x,int incx, ARRAY y,int incy);
#endif // USMV_H
