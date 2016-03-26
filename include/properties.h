#ifndef PROPERTIES_H
#define PROPERTIES_H
// **********************************************************************
//     Author : luoyulong
//
//     Date of last modification : 8.22.2012 
//      
//     Description : CONTAINS ALL CONSTANTS USED FROM THE LIBRARY
//                   CONTAINS ROUTINES FOR MANIPULATING THE ARRAYS
//                            INFOA & DESCRA OF THE DERIVED DATATYPE
//                            FOR SPARSE MATRICES
//                   get_descra:returns matrix properties stored as chars
//                   set_descra:translates integer in character descript.
//                             of sparse matrix, used by uscr
//                   get_infoa:returns matrix properties stored as ints
//                   set_infoa:sets matrix properties stored as ints
//
// **********************************************************************

#include "blas_sparse_namedconstants.h"

// *** Description of basic derived data types
#define   no_of_types  5
#define   ISP_MATRIX  0
#define   DSP_MATRIX  1
#define   SSP_MATRIX  2
#define   CSP_MATRIX  3
#define   ZSP_MATRIX  4

// *** Determine, if array indices start at 0 or 1
#define   C_BASE  0
#define   F_BASE  1
// *** Determine, if matrix is reference or copy of original data
#define   REF_OF_SOURCE  0
#define   COP_OF_SOURCE  1
// *** Determine, if the matrix is needed or its (conjugate) transpose
#define  HERMIT_MATRIX  2

#define  MEM_ALIGNMENT  64
#define  MIN_ITERATOR_NUM 4
#define  BEST_MIC_THREADS 244
#define  GPU_BLOCK_SIZE   256
#define  NUM_DIAGS_LIMIT_PARA 20
#define  FILE_NAME  "./spmv_features"


void get_descra(char* descra,int descriptor,char* message,int *ierr);
void set_descra(char* descra, int prpty,int *ierr);
void get_infoa(int* infoa,char descr, int* val,int *ierr);
void set_infoa(int* infoa, char descr,int val,int *ierr);
#endif
