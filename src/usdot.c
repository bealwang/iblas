#include "USDOT.h"
#include "comm_tools.h"
#include "malloc.h"
double BLAS_dusdot (int nz, const double* x, const int* index, const double* y, int incy, enum blas_base_type index_base)
{
  double*  zy;
  double dusdot=0;
  int i=0;

  if(nz<=0)
    dusdot=0;
  else
    {
      zy=(double*)aligned_malloc(sizeof(double)*(incy+index_base));
      //#pragma omp parallel for num_threads(dtn(M,MIN_ITERATOR_NUM))
      for(i=index_base;i<index_base+incy;i++)
        zy[i]=y[index[i]];

      for(i=index_base;i<index_base+nz;i++)
            dusdot+=x[i]*zy[i];
      aligned_free(zy);
    }
  return dusdot;
}
