#include "USSC.h"
void BLAS_dussc (int nz, const double* x, double* y, int incy, const int* indx,enum blas_base_type index_base)
{
  int i;
  if(nz>0)
    for(i=index_base;i<index_base+nz;i++)
      y[indx[i]]=x[i];
}
