#include "uscr/USCR_INSERT_ENTRIES.h"
#include "uscr/USCR_INSERT_ENTRY.h"
int BLAS_duscr_insert_entries (d_matrix *A, int nz, const double* val, const int* indx, const int* jndx)
{
  int i=0;
  int istat=-1;
  for(i=0;i<nz;i++)
    {
      istat=BLAS_duscr_insert_entry(A,val[i],indx[i],jndx[i]);
      if(istat!=0) return istat;
    }
  return istat;
}
