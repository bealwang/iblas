#ifndef USCR_BSC_H
#define USCR_BSC_H
#include "link.h"
#include "properties.h"
#include "comm_tools.h"
#include "usds.h"
dsp_linknode* duscr_bsc (int m,int n,double* val,int n_val,int* bindx,int n_bindx,int* bpntrb,int n_bpntrb ,int* bpntre,int n_bpntre,int mb,int kb,int lb,int prpty,int* istat);
#endif // USCR_BSC_H
