#ifndef USCR_INSERT_ENTRIES_H
#define USCR_INSERT_ENTRIES_H
#include "blas_enum.h"
#include "types.h"
#include "INSERTING.h"

int BLAS_iuscr_insert_entries(blas_sparse_matrix A,int nz,const int* val,const int *indx,const int *jndx);


int BLAS_suscr_insert_entries(blas_sparse_matrix A,int nz,const float* val,const int *indx,const int *jndx);


int BLAS_duscr_insert_entries(d_matrix* A,int nz,const double* val,const int* indx,const int* jndx);


int BLAS_cuscr_insert_entries(blas_sparse_matrix A,int nz,const complex_f* val,const int *indx,const int *jndx);


int BLAS_zuscr_insert_entries(blas_sparse_matrix A,int nz,const complex_d* val,const int *indx,const int *jndx);
#endif // USCR_INSERT_ENTRIES_H
