#include "mic/mic_lmbv_bsr.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'BSR'-STORAGE
//                   lmbv = Left Multiplication By Vector: y^T=x^TA
// **********************************************************************
void dlmbv_bsr (DSPMAT* mat,double* x,int n,double* y,int m,int* ierr)
{
  int base,ofs,bofs,i,pntr,mm,nn,mb,nb,nn_sq;
  char diag,type,part,store;
  int aa,xx,yy;//index for A,X,Y in Y=XA,just for easy reading code
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
  if ((mat->FIDA!=BSR_FORMAT)||(mat->M!=n)||(mat->K!=m)||(mm!=nn)) {
      *ierr = blas_error_param;
      return;
    }
  memset (y,0,sizeof(double)*m);;
  nn_sq = nn*nn;;
  if (diag=='U') { //process unstored diagonal
      if (m==n)
        memcpy(y,x,sizeof(double)*m);
      else{
          *ierr = blas_error_param;
          return;
        }
    }
  if ((type=='S')&&(!(part=='B'))&&(m==n))
    {
      if (part=='L')
        {
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  if(i==mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=i*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);


                    }
                  else if (i>mat->IA1[pntr+ofs] + ofs)
                    {
                      aa=(pntr+bofs)*nn_sq;
                      xx=i*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(pntr+bofs)*nn_sq;
                      xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                      yy=i*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  pntr++;
                }

            }
        }
      else{
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  if(i==mat->IA1[pntr+ofs] + ofs) {
                      aa=(pntr+bofs)*nn_sq;
                      xx=i*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);


                    }
                  else if (i<mat->IA1[pntr+ofs] + ofs){
                      aa=(pntr+bofs)*nn_sq;
                      xx=i*nn;
                      yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                      dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                      aa=(pntr+bofs)*nn_sq;
                      xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                      yy=i*nn;
                      dblock_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                    }
                  pntr++;
                }

            }
        }
      *ierr = 0;
    }
  else if((type=='H')&&(!(part=='B'))&&(m==n)) {
          if (part=='L')
            {
              for(i=0;i<mb;i++)
                {
                  pntr = mat->PB[i];
                  while(pntr<mat->PE[i])
                    {
                      if(i==mat->IA1[pntr+ofs] + ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=i*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);


                        }
                      else if (i>mat->IA1[pntr+ofs] + ofs)
                        {
                          aa=(pntr+bofs)*nn_sq;
                          xx=i*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=i*nn;
                          dblock_Z_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }

                }
            }
          else{
              for(i=0;i<mb;i++)
                {
                  pntr = mat->PB[i];
                  while(pntr<mat->PE[i])
                    {
                      if(i==mat->IA1[pntr+ofs] + ofs) {
                          aa=(pntr+bofs)*nn_sq;
                          xx=i*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);


                        }
                      else if (i<mat->IA1[pntr+ofs] + ofs){
                          aa=(pntr+bofs)*nn_sq;
                          xx=i*nn;
                          yy=(mat->IA1[pntr+ofs]+bofs)*nn;
                          dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);

                          aa=(pntr+bofs)*nn_sq;
                          xx=(mat->IA1[pntr+ofs]+bofs)*nn;
                          yy=i*nn;
                          dblock_Z_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
                        }
                      pntr++;
                    }

                }
            }
          *ierr = 0;
        }
      else{ //no symmetry
        // #pragma omp parallel for num_threads(dtn(mb,MIN_ITERATOR_NUM))
        //   for(i=0;i<mb;i++)
        //     {
        //       int row_begin = mat->PB[i];
        //       int row_end = mat->PE[i];
        //       int j;
        //       for(j=row_begin;j<row_end;j++)
        //         {
        //           //dblock_T_mult_vec(&mat->A[aa],&x[xx],nn,&y[yy],nn,store,ierr);
        //           int p;
        //           for (p=0;p<nn;p++){
        //             int q;
        //             for (q=0;q<nn;q++)
        //               {
        //                 if((mat->IA1[j]*nn+q < mat->K) && (i*nn+p < mat->M)){
        //                   y[mat->IA1[j]*nn+q]+=mat->A[j*nn*nn+p*nn+q]*x[i*nn+p];
        //                 }
        //               }
        //           }
        //         }
        //     }
        lmbv_bsr(mat->IA1,mat->A,mat->PB,mat->PE,mb,nb,mat->M,mat->K,mat->n_IA1,nn,x,y);
          *ierr = 0;
        }
    }// dlmbv_bsr

