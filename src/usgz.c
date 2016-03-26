#include "USGZ.h"
#include "USGA.h"
void BlAS_dusgz (int nz, double* y, int incy, double* x, const int* indx,enum blas_base_type index_base)
{
  int i;
  if(nz>0)
    {
      BLAS_dusga (nz,y,nz,x,indx,index_base);
      for(i=index_base;i<index_base+nz;i++)
        y[indx[i]]=0;
    }
}
