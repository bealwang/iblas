#ifndef CONV_TOOLS_H
#define CONV_TOOLS_H
#include "blas_sparse_namedconstants.h"
#include "types.h"
#include <malloc.h>
#include "comm_tools.h"

// pre_usconv_coo2csr
void dpre_usconv_coo2csr(int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int* PB,int* PE);

// pre_usconv_coo2csc
void dpre_usconv_coo2csc(int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int* PB,int* PE);

// pre_usconv_bco2bsc
void dpre_usconv_bco2bsc (int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int lb,int* PB,int* PE);

// pre_usconv_bco2bsr
void dpre_usconv_bco2bsr (int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int lb,int* PB,int* PE);

// pre_usconv_coo2dia
void dpre_usconv_coo2dia  (int m,int n,double** VAL,int n_VAL,int** INDX,int* JNDX,int* LDA,int* NDIAG);

// pre_usconv_dia2coo
void dpre_usconv_dia2coo (double** VAL_DIA,int** IDIAG,int n_IDIAG ,int** IA2,int LDA,int NNZ,int m);

// pre_usconv_bco2bdi
void dpre_usconv_bco2bdi (int mb,int kb,int lb,double **VAL,int** BINDX,int n_BINDX,int* BJNDX,int *BLDA,int *BNDIAG);

// pre_usconv_bdi2bco
void dpre_usconv_bdi2bco (double** VAL_DIA,int n_VAL_DIA,int** BIDIAG,int n_BIDIAG,int** IA2,int BLDA,int* BNNZ,int lb);

void dpre_usconv_coo2bco (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int row_block_size, int col_block_size,int* bnnz,int*mb,int* kb);
void dpre_usconv_bco2coo (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int* nnz,int bnnz,int lb,int mb,int kb);
void dpre_usconv_csr2bco  (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int* PB,int* PE,int row_block_size,int col_block_size,int* bnnz,int* mb,int* kb);

#endif // CONV_TOOLS_H
