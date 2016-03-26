#include "gpu/gpu_lmbv_csc.h"
#include <stdlib.h>
#include <string.h>
/*
*
*Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
*lmbv=Left Multiplication By Vector: y=xA
*/
void dlmbv_csc (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int base,ofs,j,pntr;
  char diag,type,part;
  *ierr=-1;

  if (mat->FIDA!=CSC_FORMAT||mat->M!=n||mat->K!=m)
    {
      *ierr=blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'b',&base,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return;
    }
ofs=base;
  get_descra(mat->DESCRA,'d',&diag,ierr);

  get_descra(mat->DESCRA,'t',&type,ierr);

  get_descra(mat->DESCRA,'a',&part,ierr);

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
        }
    }

  if ((type=='S'||type=='H')&&part!='B'&&m==n ) // S :symetic  B:both lower and upper are specified!B means either upper or lower are specified
    {
      if (part=='L') //lower are specified
        {
          for (j=0;j<mat->M;j++)
            {
              pntr=mat->PB[j];//begin
              while (pntr<mat->PE[j]) //end
                {
                  if (j==mat->IA1[pntr+ofs]+ofs )   // IA1 may be col idx
                    y[j]+=mat->A[pntr+ofs]*x[j];
                  else if (j<mat->IA1[pntr+ofs]+ofs ) //
                    {
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                      y[j]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                    }
                  pntr+=1;
                }
            }
        }
      else //upper
        {
          for (j=0;j<mat->M;j++)
            {
              pntr=mat->PB[j];//begin
              while (pntr<mat->PE[j]) //end
                {
                  if (j==mat->IA1[pntr+ofs]+ofs )   // IA1 may be col idx
                    y[j]+=mat->A[pntr+ofs]*x[j];//y[j]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs]
                  else if (j>mat->IA1[pntr+ofs]+ofs ) //
                    {
                      y[j]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                    }
                  pntr++;
                }
            }
        }
    }
  else //general pattern
    {
      // #pragma omp parallel for num_threads(dtn(mat->K,MIN_ITERATOR_NUM))
      // for (j=0;j<mat->K;j++)
      //   {
      //     int col_begin = mat->PB[j];
      //     int col_end = mat->PE[j];
      //     int i;
      //     double sum = 0.0;
      //     for (i=col_begin;i<col_end;i++)
      //       {
      //         sum+=mat->A[i]*x[mat->IA1[i]];
      //       }
      //     y[j] = sum;
      //   }
      lmbv_csc(mat->IA1,mat->A,mat->PB,mat->PE,mat->M,mat->K,mat->n_IA1,x,y);
      *ierr=0;
    }
}
