#ifndef USCR_CSR_H
#define USCR_CSR_H
#include "link.h"
#include "properties.h"
#include "comm_tools.h"
#include "usds.h"
dsp_linknode* duscr_csr (int m,int n,double* val,int n_val,int* indx,int n_indx,int* pntrb,int n_pntrb,int* pntre,int n_pntre,int prpty,int* istat);
#endif // USCR_CSR_H
