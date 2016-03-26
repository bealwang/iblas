#ifndef GPU_RMBV_COO_H
#define GPU_RMBV_COO_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
int drmbv_coo(DSPMAT* mat, double *x, int n,  double *y, int m, int *ierr);
#endif
