#include "mic/mic_lmbv_csr.h"
#include <stdlib.h>
#include <string.h>
/*
*
*Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
*lmbv=Left Multiplication By Vector: y=xA
*/
void dlmbv_csr(DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int base,ofs,i,pntr;
  char diag,type,part;
  *ierr=-1;
  if (mat->FIDA!=CSR_FORMAT||mat->M!=n||mat->K!=m)
    {
      *ierr=blas_error_param;
      return;
    }

  get_infoa(mat->INFOA,'b',&base,ierr);
ofs=base;
  get_descra(mat->DESCRA,'d',&diag,ierr);

  get_descra(mat->DESCRA,'t',&type,ierr);

  get_descra(mat->DESCRA,'a',&part,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return;
    }
  memset(y,0,m*sizeof(double));
  if(diag=='U') //unstored diagonal
    {
      if (m==n )
        {
          memcpy(y,x,sizeof(double)*m);
        }
      else
        {
          *ierr=blas_error_memalloc;
          return;
        }
    }
  if ((type=='S'||type=='H')&&part!='B'&&m==n ) // S :symetic  B:both lower and upper are specified!B means either upper or lower are specified
    {
      if (part=='U') //upper are specified
        {
          for (i=0;i<mat->M;i++)
            {
              pntr=mat->PB[i];//begin
              while (pntr<mat->PE[i]) //end
                {
                  if (i==mat->IA1[pntr+ofs]+ofs )   // IA1 may be col idx
                    y[i]+=mat->A[pntr+ofs]*x[i];
                  else if (i<mat->IA1[pntr+ofs]+ofs ) //
                    {
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[i];
                      y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                    }
                  pntr++;
                }
            }
        }
      else //lower
        {
          for (i=0;i<mat->M;i++)
            {
              pntr=mat->PB[i];//begin
              while (pntr<mat->PE[i]) //end
                {
                  if (i==mat->IA1[pntr+ofs]+ofs )   // IA1 may be col idx
                    y[i]+=mat->A[pntr+ofs]*x[i];
                  else if (i>mat->IA1[pntr+ofs]+ofs ) //
                    {
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[i];
                      y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                    }
                  pntr++;
                }
            }
        }
    }
  else //general pattern
    {
      // #pragma omp parallel for num_threads(dtn(mat->M,MIN_ITERATOR_NUM))
      // for (i=0;i<mat->M;i++)
      //   {
      //     int row_begin = mat->PB[i];
      //     int row_end = mat->PE[i];
      //     int j;
      //     for (j=row_begin;j<row_end;j++)
      //       {
      //         y[mat->IA1[j]]+=mat->A[j]*x[i];
      //       }
      //   }
      lmbv_csr(mat->IA1,mat->A,mat->PB,mat->PE,mat->M,mat->K,mat->n_A,x,y);

    }
}
