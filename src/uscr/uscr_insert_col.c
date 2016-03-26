#include "uscr/USCR_INSERT_COL.h"
#include "uscr/USCR_INSERT_ENTRY.h"
int BLAS_duscr_insert_col (d_matrix* A, int j, int nz, const double* val, const int* indx)
{
  int i;
  int istat=-1;

  for(i=0;i<nz;i++)
    {
      BLAS_duscr_insert_entry (A,val[i],indx[i],j);
      if(istat!=0) return istat;
    }
  return istat;
}
