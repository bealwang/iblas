#ifndef USCR_BCO_H
#define USCR_BCO_H
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'BCO'-FORMAT
// **********************************************************************
#include "link.h"
#include "properties.h"
#include "uscr/usds.h"
dsp_linknode* duscr_bco(int m,int n,double *val,int n_val,int* bindx,int n_bindx,int* bjndx,int n_bjndx,int bnnz,int mb,int kb,int lb,int prpty,int* istat);
#endif
