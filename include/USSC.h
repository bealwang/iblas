#ifndef USSC_H
#define USSC_H
// **********************************************************************
//     Author : luoyulong
//
//     Date of last modification : 8.22.2012
//
//     Description : sparse scatter
//
// **********************************************************************
#include "blas_enum.h"
//void BLAS_iussc(int nz,const ARRAY x,ARRAY y,int incy,const int *indx,enum blas_base_type);
//void BLAS_sussc(int nz,const ARRAY x,ARRAY y,int incy,const int *indx,enum blas_base_type);
void BLAS_dussc(int nz,const double* x,double* y,int incy,const int* indx,enum blas_base_type index_base);
//void BLAS_cussc(int nz,const ARRAY x,ARRAY y,int incy,const int *indx,enum blas_base_type);
//void BLAS_zussc(int nz,const ARRAY x,ARRAY y,int incy,const int *indx,enum blas_base_type);
#endif // USSC_H
