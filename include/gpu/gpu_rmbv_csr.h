#ifndef GPU_RMBV_CSR_H
#define GPU_RMBV_CSR_H
#include "types.h"
#include "properties.h"
#include "gpu_kernel_utils.h"
void drmbv_csr(DSPMAT* mat,double* x, int n, double* y, int m, int* ierr);
#endif
