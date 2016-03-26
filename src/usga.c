#include "USGA.h"

void BLAS_dusga (int nz, const double* y, int incy, double *x, const int* indx, enum blas_base_type index_base)
{
      int i;
      if(nz>0)
         for(i=index_base;i<nz+index_base;i++)
            x[i]=y[indx[i]];
}
