#ifndef GPU_LMBV_BSC_H
#define GPU_LMBV_BSC_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void dlmbv_bsc (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // LMBV_BSC_H
