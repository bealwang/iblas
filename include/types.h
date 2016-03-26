#ifndef spblas_c_types_h
#define spblas_c_types_h
#include "blas_sparse_namedconstants.h"
#include "blas_enum.h"
// **********************************************************************
//    Author : luoyulong
//
//    Date of last modification : 2012.08.22
//
//     Description : CONTAINS THE BASIC TYPES FOR SPARSE MATRICES/VECTORS
// **********************************************************************
typedef struct{
    int M;
    int N;
    int NNZ;
    int max_RD;
    double avg_RD; 
    double var_RD;   
    int Ndiags;
    double NTdiags_ratio;
    double er_dia;
    double R;
}SpFeature_CPU;

typedef struct
{
    float x;
    float y;
}complex_f;

typedef struct
{
    double x;
    double y;
}complex_d;

typedef struct
{
    int M,K;
   blas_store_format FIDA ;
    char DESCRA[11];
    int INFOA[10];
    int *A;
    int *IA1,*IA2,*PB,*PE,*BP1,*BP2;
    int n_A,n_IA1,n_IA2,n_PB,n_PE,n_BP1,n_BP2;
    int NTdiags;
    double NTdiags_ratio;
    double ER_DIA;
}ISPMAT;

typedef struct
{
    int M,K;
    blas_store_format FIDA;
    char DESCRA[11];
    int INFOA[10];
    double *A;
    int *IA1,*IA2,*PB,*PE,*BP1,*BP2;
    int n_A,n_IA1,n_IA2,n_PB,n_PE,n_BP1,n_BP2;
    SpFeature_CPU spc;
}DSPMAT;


typedef struct
{
    int M,K;
    blas_store_format FIDA ;
    char DESCRA[11];
    int INFOA[10];
    float *A;
    int *IA1,*IA2,*PB,*PE,*BP1,*BP2;
    int n_A,n_IA1,n_IA2,n_PB,n_PE,n_BP1,n_BP2;
}SSPMAT;

typedef struct
{
    int M,K;
    blas_store_format FIDA ;
    char DESCRA[11];
    int INFOA[10];
    complex_f *A;
    int *IA1,*IA2,*PB,*PE,*BP1,*BP2;
    int n_A,n_IA1,n_IA2,n_PB,n_PE,n_BP1,n_BP2;
}CSPMAT;
typedef struct
{
    int M,K;
    blas_store_format FIDA ;
    char DESCRA[11];
    int INFOA[10];
    complex_d *A;
   int *IA1,*IA2,*PB,*PE,*BP1,*BP2;
    int n_A,n_IA1,n_IA2,n_PB,n_PE,n_BP1,n_BP2;
}ZSPMAT;

#endif
