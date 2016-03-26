#ifndef GPU_RMBV_DIA_H
#define GPU_RMBV_DIA_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void drmbv_dia (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr);
#endif // RMBV_DIA_H
