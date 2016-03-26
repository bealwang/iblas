#ifndef GPU_LMBV_BSR_H
#define GPU_LMBV_BSR_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void dlmbv_bsr (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // LMBV_BSR_H
