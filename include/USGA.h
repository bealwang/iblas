#ifndef USGA_H
#define USGA_H
#include "blas_enum.h"
// **********************************************************************
//     Author : luoyulong
//
//     Date of last modification : 8.22.2012
//
//     Description : sparse gather
//
// **********************************************************************
//void BLAS_iusga(int nz,const ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
//void BLAS_susga(int nz,const ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
void BLAS_dusga(int nz,const double* y,int incy,double* x,const int* indx,enum blas_base_type index_base);
//void BLAS_cusga(int nz,const ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
//void BLAS_zusga(int nz,const ARRAY y,int incy,ARRAY x,const int *indx,enum blas_base_type index_base);
#endif // USGA_H
