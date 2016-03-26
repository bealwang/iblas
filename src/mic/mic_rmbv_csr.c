#include "mic/mic_rmbv_csr.h"
#include <stdlib.h>
#include <string.h>
/*
*
*Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
*rmbv=Right Multiplication By Vector: y=Ax
*/
void drmbv_csr (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int base,ofs,i,pntr;
  char diag,type,part;
  *ierr=-1;
  
  if (mat->FIDA!=CSR_FORMAT||mat->M!=m||mat->K!=n)
    {
      *ierr=blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'b',&base,ierr);

  get_descra(mat->DESCRA,'d',&diag,ierr);

  get_descra(mat->DESCRA,'t',&type,ierr);

  get_descra(mat->DESCRA,'a',&part,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return;
    }

  ofs=1-base;
  ofs=0;


  memset(y,0,m*sizeof(double));
  if(diag=='U') //unstored diagonal
    {
      if (m==n)
        {
          memcpy(y,x,sizeof(double)*m);
        }
      else
        {
          *ierr=blas_error_memalloc;
          return;
        }
    }
  if ((type=='S'||type=='H')&&part!='B'&&m==n) // S :symetic  B:both lower and upper are specified!B means either upper or lower are specified
    {
      if (part=='U') //upper are specified
        {
          for (i=0;i<mat->M;i++)
            {
              pntr=mat->PB[i]; //begin
              while (pntr<mat->PE[i]) //end
                {
                  if (i==mat->IA1[pntr+ofs]+ofs)   // IA1 may be col idx
                    y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                  else if (i<mat->IA1[pntr+ofs]+ofs)
                    {
                      y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[i];
                    }
                  pntr++;
                }
            }
          *ierr=0;
        }
      else //lower
        {
          for (i=0;i<mat->M;i++)
            {
              pntr=mat->PB[i]; //begin
              while (pntr<mat->PE[i]) //end
                {
                  if (i==mat->IA1[pntr+ofs]+ofs)   // IA1 may be col idx
                    y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                  else if (i>mat->IA1[pntr+ofs]+ofs)
                    {
                      y[i]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                      y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[i];
                    }
                  pntr++;
                }
            }
          *ierr=0;
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
      //     double sum = 0.0;
      //     for(j=row_begin;j<row_end;j++)
      //       {
      //         sum+=mat->A[j]*x[mat->IA1[j]];
      //       }
      //     y[i] = sum;
      //   }
      rmbv_csr(mat->IA1,mat->A,mat->PB,mat->PE,mat->M,mat->K,mat->n_A,x,y);
      *ierr=0;
    }
}
