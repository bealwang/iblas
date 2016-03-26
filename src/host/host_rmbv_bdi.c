#include "host/host_rmbv_bdi.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'BDI'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
// **********************************************************************
void drmbv_bdi (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int i,j,mm,nn,nn_sq,mb,nb;
  int blda,nbdiag,start_a,end_a,start_x,start_y;
  int aa,xx,yy;
  char diag,type,part,store;
  *ierr = -1;


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
  get_infoa(mat->INFOA,'f',&blda,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'g',&nbdiag,ierr);
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

  mb = (mat->M+nn-1)/mm;
  nb = (mat->K+nn-1)/nn;

  if ((mat->FIDA!=BDI_FORMAT)||(mat->M!=m)||(mat->K!=n)||(mm!=nn)) {
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
          for(i=0;i<nbdiag;i++)
            {// in this situation,start_x=0
              start_x =(mat->IA1[i])>0?(mat->IA1[i]):0;
              start_y = (-mat->IA1[i])>0?(-mat->IA1[i]):0;
              if (mat->IA1[i]>0)
                {
                  start_a=i*blda;
                  //here end_a equal the length of diagonal
                  end_a=(blda-mat->IA1[i])+((mat->K/nn)-blda);

                  for(j=0;j<end_a;j++)
                    {

                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*blda;
                  end_a=blda;

                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
          for(i=0;i<nbdiag;i++)
            {// in this situation,start_y=0
              start_x =(mat->IA1[i])>0?(mat->IA1[i]):0;
              start_y = (-mat->IA1[i])>0?(-mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(blda+mat->IA1[i])+((mat->M/nn)-blda);
                  start_a=(i+1)*blda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*blda;
                  end_a=blda;

                  for(j=0;j<blda;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
          for(i=0;i<nbdiag;i++)
            {// in this situation,start_x=0
              start_x =(mat->IA1[i])>0?(mat->IA1[i]):0;
              start_y = (-mat->IA1[i])>0?(-mat->IA1[i]):0;
              if (mat->IA1[i]>0)
                {
                  start_a=i*blda;
                  //here end_a equal the length of diagonal
                  end_a=(blda-mat->IA1[i])+((mat->K/nn)-blda);

                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*blda;
                  end_a=blda;

                  for(j=0;j<blda;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                }
              else
                {
                  continue;
                }
            }
        }
      else//lower
        {
          for(i=0;i<nbdiag;i++)
            {// in this situation,start_y=0
              start_x =(mat->IA1[i])>0?(mat->IA1[i]):0;
              start_y = (-mat->IA1[i])>0?(-mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(blda+mat->IA1[i])+((mat->M/nn)-blda);
                  start_a=(i+1)*blda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_H_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                    }
                }
              else if (mat->IA1[i]==0)
                {
                  start_a=i*blda;
                  end_a=blda;

                  for(j=0;j<blda;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

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
      #pragma omp parallel for num_threads(dtn(mb,MIN_ITERATOR_NUM))
      for (i=0; i<mb; ++i) {
        int j;
        for (j=0; j<nbdiag; ++j){
          int offset = mat->IA1[j];
          int istart = offset >= 0?0:(-offset);
          int jstart = offset >= 0?offset:0;
          int N = 0; int v_offset = 0; int J_offset = 0;

          if((mb-istart) >= (nb-jstart)){
            N = nb - jstart;
            J_offset = jstart - istart;
            v_offset = j*blda - istart;
          }else{
            N = mb - istart;
            J_offset = jstart - istart;
            v_offset = j*blda + blda - mb;
          }

          int iend = istart + N;
          
          if ((i >= istart) && (i < iend)){
            int p;
            for (p=0;p<mm;p++){
              int q;
              for (q=0;q<nn;q++){
                if((i*mm+p < mat->M) && ((i+J_offset)*nn+q < mat->K)){
                  y[i*mm+p] += mat->A[(i+v_offset)*mm*nn+nn*p+q] * x[(i+J_offset)*nn+q];
                }
              }
            }
          }
          if(i >= iend){
            continue;
          }
        }
      }




      // for(i=0;i<nbdiag;i++)
      //   {
      //     start_x=(mat->IA1[i])>0?(mat->IA1[i]):0;
      //     start_y=(-mat->IA1[i])>0?(-mat->IA1[i]):0;

      //     if (mat->IA1[i]>(mat->K/nn)-blda){
      //         start_a=i*blda;
      //         end_a=(blda-mat->IA1[i])+((mat->K/nn)-blda);
      //       }
      //     else if (mat->IA1[i]<-(mat->M/nn)+blda)
      //       {
      //         end_a=(blda+mat->IA1[i])+((mat->M/nn)-blda);
      //         start_a=(i+1)*blda-end_a;
      //       }
      //     else
      //       {
      //         start_a=i*blda;
      //         end_a=blda;
      //       }
      //     for(j=0;j<end_a;j++)
      //       {
      //         aa=(start_a+j)*nn_sq;
      //         xx=(start_x+j)*nn;
      //         yy=(start_y+j)*nn;
      //         dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
      //       }


      //   }
      *ierr = 0;
    }

}// drmbv_bdi

