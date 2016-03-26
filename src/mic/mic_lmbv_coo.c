#include "mic/mic_lmbv_coo.h"
#include <stdio.h>
#include <string.h>
/*
 *
 *Description : PERFORMS MV MULT. WITH MATRIX IN 'CSR'-STORAGE
 *rmbv=Right Multiplication By Vector: y=Ax
 */
int dlmbv_coo(DSPMAT* mat,double* x,int n,double* y,int m,int *ierr)
{
  int nnz,base,ofs,i;
  int part;
  char diag,type;
  *ierr=-1;
  if (mat->FIDA != COO_FORMAT||mat->M!=n||mat->K!=m)
    {
      *ierr=blas_error_param;
      return 0;
    }
  get_infoa(mat->INFOA,'b',&base,ierr);
  ofs=base;

  get_infoa(mat->INFOA,'n',&nnz,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return 0;
    }
  get_descra (mat->DESCRA,'d',&diag,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return 0;
    }
  get_descra (mat->DESCRA,'t',&type,ierr);
  if (*ierr!=0)
    {
      *ierr=blas_error_param;
      return 0;
    }
  get_infoa(mat->INFOA,'a',&part,ierr);
  memset(y,0,m*sizeof(double));//double array only can be set to zero

  if(diag=='U') //unstored diagonal
    {
      if (m==n )
        {
          memcpy(y,x,sizeof(double)*m);
        }
      else
        {
          *ierr=blas_error_param;
          return 0;
          //error
        }
    }
  if (type=='S' && part !='B' && m==n ) // S :symetic  B:both lower and upper are specified  !B means either upper or lower are specified
    {
      if (part=='L') //upper are specified
        {
          for (i=0; i < nnz; i++)
            {
              if (mat->IA1[i]==mat->IA2[i])
                y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs];
              else if (mat->IA1[i] > mat->IA2[i])
                {
                  y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
                  y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs]; //lower part

                }
            }
        }
      else //lower are specified
        {
          for (i=0; i < nnz; i++)
            {
              if (mat->IA1[i]==mat->IA2[i])
                y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
              else if (mat->IA1[i] < mat->IA2[i])
                {
                  y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
                  y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs]; //lower part

                }
            }

        }

      //
    }
  else if (type=='H' && part !='B' && m==n ) //Hermitian matrix ,whose tranpa equal to itself
    {
      if (part=='L') //upper are specified
        {
          for (i=0; i < nnz; i++)
            {
              if (mat->IA1[i]==mat->IA2[i])
                y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
              else if (mat->IA1[i] > mat->IA2[i])
                {
                  y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
                  y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs]; //lower part

                }
            }
        }
      else //lower are specified
        {
          for (i=0; i < nnz; i++)
            {
              if (mat->IA1[i]==mat->IA2[i])
                y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
              else if (mat->IA1[i]<mat->IA2[i])
                {
                  y[mat->IA1[i]+ofs]+=mat->A[i]*x[mat->IA2[i]+ofs];
                  y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs]; //lower part

                }
            }

        }
    }
  else //general pattern
    {
      // for (i=0; i < nnz; i++)
      //   {
      //     y[mat->IA2[i]+ofs]+=mat->A[i]*x[mat->IA1[i]+ofs];
      //   }
      lmbv_coo(mat->IA1,mat->IA2,mat->A,mat->M,mat->K,nnz,x,y);
      *ierr=0;
    }
    return 1;
}
