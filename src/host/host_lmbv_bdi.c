#include "host/host_lmbv_bdi.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'BDI'-STORAGE
//                   lmbv = Left Multiplication By Vector: y^T=x^TA
// **********************************************************************
void dlmbv_bdi (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{

  int i,j,mm,nn,nn_sq,mb,nb;
  int blda,nbdiag,start_a,end_a,start_x,start_y;
  char diag,type,part,store;
  int aa,xx,yy;
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

  if ((mat->FIDA!=BDI_FORMAT)||(mat->M!=n)||((mat->K)!=m)||(mm!=nn)) {
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
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(blda+mat->IA1[i])+((mat->M/nn)-blda);
                  start_a=(i+1)*blda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_Z_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
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
              start_x =(-mat->IA1[i])>0?(-mat->IA1[i]):0;
              start_y = (mat->IA1[i])>0?(mat->IA1[i]):0;
              if (mat->IA1[i]<0) {

                  end_a=(blda+mat->IA1[i])+((mat->M/nn)-blda);
                  start_a=(i+1)*blda-end_a;
                  for(j=0;j<end_a;j++)
                    {
                      aa=(start_a+j)*nn_sq;
                      xx=(start_x+j)*nn;
                      yy=(start_y+j)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(start_a+j)*nn_sq;
                      xx=(start_y+j)*nn;
                      yy=(start_x+j)*nn;
                      dblock_Z_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

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
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

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

      #pragma omp parallel for num_threads(dtn(nb,MIN_ITERATOR_NUM))
      for (i=0; i<nb; ++i) {
        int j;
        for (j=0; j<nbdiag; ++j){
          int offset = mat->IA1[j];
          int istart = offset >= 0?0:(-offset);
          int jstart = offset >= 0?offset:0;
          int N = 0; int v_offset = 0; int I_offset = 0;

          if((mb-istart) > (nb-jstart)){
            N = nb - jstart;
            I_offset = istart - jstart;
            v_offset = j*blda - jstart;
          }else{
            N = mb - istart;
            I_offset = istart - jstart;
            v_offset = j*blda + blda - jstart - N;
          }

          int jend = jstart + N;
          
          if ((i >= jstart) && (i < jend)){
            int p;
            for (p=0;p<nn;p++){
              int q;
              for (q=0;q<nn;q++){
                if((i*nn+q < mat->K) && ((i+I_offset)*mm+p < mat->M)){
                  y[i*nn+q] += mat->A[(i+v_offset)*mm*nn+nn*p+q] * x[(i+I_offset)*mm+p];
                }
              }
            }
          }
          if(i >= jend){
            continue;
          }
        }
      }



      // for(i=0;i<nbdiag;i++)
      //   {
      //     int begin_a,end,k;
      //     int begin_x=(-mat->IA1[i])>0?(-mat->IA1[i]):0;
      //     int begin_y=(mat->IA1[i])>0?(mat->IA1[i]):0;
      //     if (mat->IA1[i]>((mat->K+nn-1)/nn)-blda){
      //         begin_a=i*blda;
      //         end=(blda-mat->IA1[i])+(((mat->K+nn-1)/nn)-blda);
      //       }
      //     else if (mat->IA1[i]<-((mat->M+nn-1)/nn)+blda)
      //       {
      //         end=(blda+mat->IA1[i])+(((mat->M+nn-1)/nn)-blda);
      //         begin_a=(i+1)*blda-end;
      //       }
      //     else
      //       {
      //         begin_a=i*blda;
      //         end=blda;
      //       }
      //     for(k=0;k<end;k++){
      //         //dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
      //         int p;
      //         for (p=0;p<nn;p++){
      //           int q;
      //           for (q=0;q<nn;q++){
      //               y[(begin_y+k)*nn+q]+=mat->A[(begin_a+k)*nn*nn+p*nn+q]*x[(begin_x+k)*nn+p];
      //           }
      //         }
      //     }
      //   }
      *ierr = 0;
    }

}// dlmbv_bdi

