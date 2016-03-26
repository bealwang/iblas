#include "gpu/gpu_rmbv_bsc.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSC'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
// **********************************************************************
void drmbv_bsc (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{

  int base,ofs,bofs,j,pntr,mm,nn,mb,nb,nn_sq;
  char diag,type,part,store;
  int aa,xx,yy;
  *ierr = -1;


  get_infoa(mat->INFOA,'b',&base,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  ofs=base;
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
  get_infoa(mat->INFOA,'f',&mb,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'g',&nb,ierr);
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
  if ((mat->FIDA!=BSC_FORMAT)||(mat->M!=m)||(mat->K!=n)||(mm!=nn))
    {
      *ierr = blas_error_param;
      return;
    }
  memset (y,0,sizeof(double)*m);
  nn_sq = nn*nn;
  if (diag=='U') { //process unstored diagonal
      if (m==n)
        memcpy (y,x,sizeof(double)*m);
      else
        {
          *ierr = blas_error_param;
          return;
        }
    }
  if ((type=='S')&&(!(part=='B'))&&(m==n))
    {
      if (part=='U')//upper
        {
          for(j=0;j<nb;j++)
            {
              pntr = mat->PB[j];
              while(pntr<mat->PE[j])
                {
                  if(j==mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=j*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  else
                    {
                      if (j>mat->IA1[pntr+ofs]+ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=j*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=j*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }
                }
            }
        }
      else //lower
        {
          for(j=0;j<nb;j++)
            {
              pntr = mat->PB[j];
              while(pntr<mat->PE[j])
                {
                  if(j==mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=j*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  else
                    {
                      if (j<mat->IA1[pntr+ofs]+ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=j*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=j*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }
                }
            }
        }
      *ierr=0;
    }
  else  if((type=='H')&&(!(part=='B'))&&(m==n))
    {
      if (part=='U')//upper
        {
          for(j=0;j<nb;j++)
            {
              pntr = mat->PB[j];
              while(pntr<mat->PE[j])
                {
                  if(j==mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=j*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  else
                    {
                      if (j>mat->IA1[pntr+ofs]+ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=j*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=j*nn;
                          dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }
                }
            }
        }
      else //lower
        {
          for(j=0;j<nb;j++)
            {
              pntr = mat->PB[j];
              while(pntr<mat->PE[j])
                {
                  if(j==mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=j*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  else
                    {
                      if (j<mat->IA1[pntr+ofs]+ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=j*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=j*nn;
                          dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }
                }
            }

        }
      *ierr=0;
    }
  else
    { //no symmetry
      #pragma omp parallel for num_threads(dtn(nb,MIN_ITERATOR_NUM))
      for(j=0;j<nb;j++)
        { 
          int col_begin = mat->PB[j];
          int col_end = mat->PE[j];
          int i;
          for(i=col_begin;i<col_end;i++)
            {
              int p;
              //dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
              for (p=0;p<nn;p++){
                int q;
                for(q=0;q<nn;q++){
                  if((j*nn+q < mat->K) && (mat->IA1[i]*nn+p < mat->M)){
                    y[mat->IA1[i]*nn+p]+=mat->A[i*nn*nn+p*nn+q]*x[j*nn+q];
                  }
                }
              }
            }
        }
      // rmbv_bsc(mat->IA1,mat->A,mat->PB,mat->PE,mb,nb,mat->M,mat->K,mat->n_IA1,nn,x,y);
      *ierr = 0;
    }
}// drmbv_bsc
