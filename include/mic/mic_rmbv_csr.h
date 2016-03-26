#ifndef MIC_RMBV_CSR_H
#define MIC_RMBV_CSR_H
#include "types.h"
#include "properties.h"
#include "mic_kernel_utils.h"
void drmbv_csr(DSPMAT* mat,double* x, int n, double* y, int m, int* ierr);
#endif
