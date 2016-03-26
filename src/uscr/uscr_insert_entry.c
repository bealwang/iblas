#include "uscr/USCR_INSERT_ENTRY.h"
int BLAS_duscr_insert_entry  (d_matrix* A,double val,int i,int j)
{
d_matrix* pmatrix;
int istat=-1;
 pmatrix=daccess_matrix(A);
switch (pmatrix->format)
{
  case 'b':
    dINS_bl_entr (pmatrix,val,i,j,&istat);
  case'v':
     dINS_varbl_entr(pmatrix,val,i,j,&istat);
  case 'n':
    dINS_entry (pmatrix,val,i,j,&istat);
  default:
   	istat=blas_error_param;
  	return istat;
}
return istat;
}//end subroutine duscr_insert_entry
