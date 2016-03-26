#ifndef GPU_LMBV_DIA_H
#define GPU_LMBV_DIA_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void dlmbv_dia (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // LMBV_DIA_H
