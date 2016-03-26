// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'COO'-FORMAT
// **********************************************************************
#include "link.h"
#include "properties.h"
#include "comm_tools.h"
#include "usds.h"
// **********************************************************************
// **********************************************************************
dsp_linknode* duscr_coo(int m,int n,double* val,int n_val,int* indx,int* jndx,int nnz,int prpty,int* istat);
