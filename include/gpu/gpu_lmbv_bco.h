#ifndef GPU_LMBV_BCO_H
#define GPU_LMBV_BCO_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void dlmbv_bco (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // LMBV_BCO_H
