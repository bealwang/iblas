#ifndef USGZ_H
#define USGZ_H
// **********************************************************************
//     Author : luoyulong
//
//     Date of last modification : 8.22.2012
//
//     Description : sparse gather and zero
//
// **********************************************************************
#include "blas_enum.h"
//void BlAS_iusgz(int nz,ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
//void BlAS_susgz(int nz,ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
void BlAS_dusgz(int nz,double* y,int incy,double* x,const int* indx,enum blas_base_type index_base);
//void BlAS_cusgz(int nz,ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
//void BlAS_zusgz(int nz,ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
#endif // USGZ_H
