#include "USAXPY.h"

void BLAS_dusaxpy (int nz, double alpha, const double* x, const int* indx, double* y, int incy, enum blas_base_type index_base)
{
  int i;
  if(nz>0)
      for(i=index_base;i<index_base+nz;i++)
        y[indx[i]]+=x[i]*alpha;
}

