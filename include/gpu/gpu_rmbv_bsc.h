#ifndef GPU_RMBV_BSC_H
#define GPU_RMBV_BSC_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"

void drmbv_bsc (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // RMBV_BSC_H
