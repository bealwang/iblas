#ifndef USAXPY_H
#define USAXPY_H
// **********************************************************************
//     Author : luoyulong
//
//     Date of last modification : 8.22.2012
//
//     Description : sparse vector update
//
// **********************************************************************
#include "blas_enum.h"

//void BLAS_iusaxpy(int nz,SCALAR_IN alpha,const ARRAY ,const int *indx,ARRAY y,int incy,enum blas_base_type index_base);

//void BLAS_susaxpy(int nz,SCALAR_IN alpha,const ARRAY ,const int *indx,ARRAY y,int incy,enum blas_base_type index_base);

void BLAS_dusaxpy (int nz, double alpha, const double* x, const int* indx, double* y, int incy,  enum blas_base_type  index_base);

//void BLAS_cusaxpy(int nz,SCALAR_IN alpha,const ARRAY ,const int *indx,ARRAY y,int incy,enum blas_base_type index_base);

//void BLAS_zusaxpy(int nz,SCALAR_IN alpha,const ARRAY ,const int *indx,ARRAY y,int incy,enum blas_base_type index_base);
#endif // USAXPY_H
