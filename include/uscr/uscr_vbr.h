#ifndef USCR_VBR_H
#define USCR_VBR_H
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : CREATION ROUTINE FOR MATRIX HANDLE FROM 'VBR'-FORMAT
// **********************************************************************
#include "link.h"
#include "properties.h"
#include "usds.h"
#include "comm_tools.h"
// **********************************************************************
// **********************************************************************
//the size of index is block number+1!
dsp_linknode* duscr_vbr(int m,int n,double* val,int n_val,int* indx,int n_indx,int* bindx,int n_bindx,int* rpntr,int n_rpntr,int* cpntr,int n_cpntr,int* bpntrb,int n_bpntrb,int* bpntre,int n_bpntre,int mb,int kb,int prpty,int* istat);
#endif
