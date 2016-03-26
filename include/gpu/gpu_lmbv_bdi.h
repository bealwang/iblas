#ifndef GPU_LMBV_BDI_H
#define GPU_LMBV_BDI_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void dlmbv_bdi (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // LMBV_BDI_H
