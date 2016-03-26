#ifndef GPU_RMBV_BSR_H
#define GPU_RMBV_BSR_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void drmbv_bsr (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // RMBV_BSR_H
