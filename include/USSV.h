#ifndef USSM_H
#define USSM_H
#include "blas_enum.h"
#include "../types.h"

int BLAS_iussv(enum blas_trans_type transt,SCALAR_IN alpha,blas_sparse_matrix T,ARRAY x,int incx);
int BLAS_sussv(enum blas_trans_type transt,SCALAR_IN alpha,blas_sparse_matrix T,ARRAY x,int incx);
int BLAS_dussv(enum blas_trans_type transt,SCALAR_IN alpha,blas_sparse_matrix T,ARRAY x,int incx);
int BLAS_cussv(enum blas_trans_type transt,SCALAR_IN alpha,blas_sparse_matrix T,ARRAY x,int incx);
int BLAS_zussv(enum blas_trans_type transt,SCALAR_IN alpha,blas_sparse_matrix T,ARRAY x,int incx);


#endif // USSM_H

