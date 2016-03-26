#include "mic/mic_rmbv_dia.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'DIA'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
// **********************************************************************
void drmbv_dia (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int i,j;
  int lda,ndiag,start_a,end_a,start_x,start_y;
  char diag,type,part;
  *ierr = -1;

  if ((mat->FIDA!=DIA_FORMAT)||(mat->M!=m)||(mat->K!=n)) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'d',&lda,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'e',&ndiag,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'d',&diag,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'t',&type,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'a',&part,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  memset (y,0,sizeof(double)*m);
  if (diag=='U') { //process unstored diagonal
      if (m==n)
        memcpy (y,x,sizeof(double)*m);
      else{
          *ierr = blas_error_param;
          return;
        }
    }
  if ((type=='S')&&(!(part=='B'))&&(m==n)) {
      if (part=='U') {
          for(i=0;i<ndiag;i++)
            {// in this situation,start_x=0
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]>0)
                {
                  start_a=i*lda;
                  //here end_a equal the length of diagonal
                  end_a=(lda-mat->IA1[i])+(mat->K-lda);

                  for(j=0;j<end_a;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                      y[start_x+j]+=mat->A[start_a+j]*x[start_y+j];
                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*lda;
                  end_a=lda;

                  for(j=0;j<lda;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                    }
                }
              else
                {
                  continue;
                }
            }
        }
      else
        {
          for(i=0;i<ndiag;i++)
            {// in this situation,start_y=0
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(lda+mat->IA1[i])+(mat->M-lda);
                  start_a=(i+1)*lda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                      y[start_x+j]+=mat->A[start_a+j]*x[start_y+j];

                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*lda;
                  end_a=lda;

                  for(j=0;j<lda;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                    }
                }
              else
                {
                  continue;
                }
            }
        }
      *ierr = 0;
    }
  else if ((type=='H')&&(!(part=='B'))&&(m==n))
    {
      if (part=='U') {
          for(i=0;i<ndiag;i++)
            {// in this situation,start_x=0
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]>0)
                {
                  start_a=i*lda;
                  //here end_a equal the length of diagonal
                  end_a=(lda-mat->IA1[i])+(mat->K-lda);

                  for(j=0;j<end_a;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                      y[start_x+j]+=mat->A[start_a+j]*x[start_y+j];
                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*lda;
                  end_a=lda;

                  for(j=0;j<lda;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                    }
                }
              else
                {
                  continue;
                }
            }
        }
      else
        {
          for(i=0;i<ndiag;i++)
            {// in this situation,start_y=0
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(lda+mat->IA1[i])+(mat->M-lda);
                  start_a=(i+1)*lda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                      y[start_x+j]+=mat->A[start_a+j]*x[start_y+j];

                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*lda;
                  end_a=lda;

                  for(j=0;j<lda;j++)
                    {
                      y[start_y+j]+=mat->A[start_a+j]*x[start_x+j];
                    }
                }
              else
                {
                  continue;
                }
            }
        }
      *ierr = 0;
    }
  else
    { //no symmetry
      rmbv_dia(mat->IA1,mat->A,mat->M,mat->K,ndiag,lda,x,y);
      *ierr = 0;
    }

}// drmbv_dia

