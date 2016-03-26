#ifndef USSM_H
#define USSM_H
#include "blas_enum.h"
int BLAS_iussm(enum blas_order_type order,enum blas_trans_type transt,int nrhs,SCALAR_IN alpha,blas_sparse_matrix T, ARRAY B,int ldb);

int BLAS_sussm(enum blas_order_type order,enum blas_trans_type transt,int nrhs,SCALAR_IN alpha,blas_sparse_matrix T, ARRAY B,int ldb);

int BLAS_dussm(enum blas_order_type order,enum blas_trans_type transt,int nrhs,SCALAR_IN alpha,blas_sparse_matrix T, ARRAY B,int ldb);

int BLAS_cussm(enum blas_order_type order,enum blas_trans_type transt,int nrhs,SCALAR_IN alpha,blas_sparse_matrix T, ARRAY B,int ldb);

int BLAS_zussm(enum blas_order_type order,enum blas_trans_type transt,int nrhs,SCALAR_IN alpha,blas_sparse_matrix T, ARRAY B,int ldb);
#endif // USSM_H
