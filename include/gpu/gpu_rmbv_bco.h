#ifndef GPU_RMBV_BCO_H
#define GPU_RMBV_BCO_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void drmbv_bco (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // RMBV_BCO_H
