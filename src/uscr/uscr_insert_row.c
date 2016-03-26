#include "uscr/USCR_INSERT_ROW.h"
#include "uscr/USCR_INSERT_ENTRY.h"
int BLAS_duscr_insert_row (d_matrix* A, int i, int nz, const double *val, const int *jndx)
{
  int k=0;
  int istat=-1;
  for(k=0;k<nz;k++)
    {
      istat=BLAS_duscr_insert_entry(A,val[k],i,jndx[k]);
      if(istat!=0) return istat;
    }
  return istat;
}
