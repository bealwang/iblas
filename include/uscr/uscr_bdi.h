#ifndef USCR_BDI_H
#define USCR_BDI_H
#include "link.h"
#include "properties.h"
#include "comm_tools.h"
#include "usds.h"
dsp_linknode* duscr_bdi (int m,int n,double* val,int n_val,int blda,int* ibdiag,int n_ibdiag,int nbdiag,int mb,int kb,int lb,int prpty,int* istat);
#endif // USCR_BDI_H
