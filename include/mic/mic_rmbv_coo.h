#ifndef MIC_RMBV_COO_H
#define MIC_RMBV_COO_H
#include "types.h"
#include "properties.h"
#include "mic_kernel_utils.h"

int drmbv_coo(DSPMAT* mat, double *x, int n,  double *y, int m, int *ierr);
#endif
