#include "host/host_rmbv_bco.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'BCO'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
// **********************************************************************
void drmbv_bco (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int base,bofs,i,mm,nn,nnz,nn_sq;
  char diag,type,part,store;
  int aa,xx,yy;
  *ierr = -1;

  get_infoa(mat->INFOA,'b',&base,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  bofs=-base;
  get_infoa(mat->INFOA,'d',&mm,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'e',&nn,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'n',&nnz,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'d',&diag,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'f',&store,ierr);
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
  if ((mat->FIDA!=BCO_FORMAT)||(mat->M!=m)||(mat->K!=n)||(mm!=nn)) {
      *ierr = blas_error_param;
      return;
    }
  memset (y,0,sizeof(double)*m);
  nn_sq = nn*nn;
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
          for(i=0;i<nnz;i++)
            {
              if(mat->IA1[i]==mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
              else if (mat->IA1[i]<mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                  aa=i*nn_sq;
                  xx=(mat->IA1[i]+bofs)*nn;
                  yy=(mat->IA2[i]+bofs)*nn;
                  dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
            }
        }
      else
        {
          for(i=0;i<nnz;i++)
            {
              if(mat->IA1[i]==mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
              else if (mat->IA1[i]<mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                  aa=i*nn_sq;
                  xx=(mat->IA1[i]+bofs)*nn;
                  yy=(mat->IA2[i]+bofs)*nn;
                  dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
            }
        }
      *ierr = 0;
    }
  else if ((type=='H')&&(!(part=='B'))&&(m==n))
    {
      if (part=='U') {//upper
          for(i=0;i<nnz;i++)
            {
              if(mat->IA1[i]==mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
              else if (mat->IA1[i]<mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                  aa=i*nn_sq;
                  xx=(mat->IA1[i]+bofs)*nn;
                  yy=(mat->IA2[i]+bofs)*nn;
                  dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
            }
        }
      else//lower
        {
          for(i=0;i<nnz;i++)
            {
              if(mat->IA1[i]==mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
              else if (mat->IA1[i]<mat->IA2[i])
                {
                  aa=i*nn_sq;
                  xx=(mat->IA2[i]+bofs)*nn;
                  yy=(mat->IA1[i]+bofs)*nn;
                  dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                  aa=i*nn_sq;
                  xx=(mat->IA1[i]+bofs)*nn;
                  yy=(mat->IA2[i]+bofs)*nn;
                  dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                }
            }
        }
      *ierr = 0;
    }
  else
    { //no symmetry
      #pragma omp parallel for num_threads(dtn(nnz,MIN_ITERATOR_NUM))
      for(i=0;i<nnz;i++)
        {
          int p;
          //dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
          for (p=0;p<nn;p++){
            int q;
            for(q=0;q<nn;q++){
              if((mat->IA1[i]*nn+p < mat->M) && (mat->IA2[i]*nn+q < mat->K)){
                y[mat->IA1[i]*nn+p]+=mat->A[i*nn*nn+p*nn+q]*x[mat->IA2[i]*nn+q];
              }
            }
          }
        }
      *ierr = 0;
    }
}// drmbv_bco
