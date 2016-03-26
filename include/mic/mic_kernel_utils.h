#ifndef MIC_KERNEL_UTILS_H
#define MIC_KERNEL_UTILS_H
#include "properties.h"
void lmbv_coo(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y);
void rmbv_coo(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y);
void lmbv_csr(int* col,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y);
void rmbv_csr(int* col,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y);
void lmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y);
void rmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y);
void lmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y);
void rmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y);

void lmbv_bco(int* row,int* col,double* val,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void rmbv_bco(int* row,int* col,double* val,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void lmbv_bsr(int* col,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void rmbv_bsr(int* col,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void lmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void rmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y);
void lmbv_bdi(int* bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y);
void rmbv_bdi(int* bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y);
#endif