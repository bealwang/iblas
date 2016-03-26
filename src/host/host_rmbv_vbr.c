#include "host/host_rmbv_vbr.h"
#include "dense.h"
#include <string.h>
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'VBR'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
// **********************************************************************
void drmbv_vbr (DSPMAT* mat,double* x,int n,double *y,int m,int *ierr)
{


  int base,ofs,i,pntr,mb,nb;
  int start_a,start_x, len_x,start_y,len_y;
  char diag,type,part,store;
  *ierr = -1;


  if ((mat->FIDA!=VBR_FORMAT)||(mat->M!=m)||(mat->K!=n)) {
      *ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'b',&base,ierr);
  if (*ierr!=0) {
      *ierr = blas_error_param;
      return;
    }
  ofs=base;
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
  memset (y,0,sizeof(double)*m);
  //      start_a = -1
  //      end_a = -1
  //      start_x = -1
  //      end_x = -1
  //      start_y = -1
  //      end_y = -1
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
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  start_y=mat->BP1[i]+ofs;
                  len_y=mat->BP1[i+1]-mat->BP1[i];

                  //IA1 block colume id to block id
                  start_x = mat->BP2[mat->IA1[pntr+ofs]+ofs] + ofs;
                  len_x = mat->BP2[mat->IA1[pntr+ofs]+ofs+1]-mat->BP2[mat->IA1[pntr+ofs]+ofs];

                  start_a = mat->IA2[pntr+ofs] + ofs;


                  if(i==mat->IA1[pntr+ofs] + ofs)
                    dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                  else  if (i<mat->IA1[pntr+ofs] + ofs)
                    {
                      dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                      dblock_T_mult_vec (&(mat->A[start_a]),&(x[start_y]),len_y,&(y[start_x]),len_x,store,ierr);
                    }
                  pntr++;
                }
            }
        }else{
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  start_y=mat->BP1[i]+ofs;
                  len_y=mat->BP1[i+1]-mat->BP1[i];

                  //IA1 block colume id to block id
                  start_x = mat->BP2[mat->IA1[pntr+ofs]+ofs] + ofs;
                  len_x = mat->BP2[mat->IA1[pntr+ofs]+ofs+1]-mat->BP2[mat->IA1[pntr+ofs]+ofs];

                  start_a = mat->IA2[pntr+ofs] + ofs;


                  if(i==mat->IA1[pntr+ofs] + ofs)
                    dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                  else  if (i>mat->IA1[pntr+ofs] + ofs)
                    {
                      dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                      dblock_T_mult_vec (&(mat->A[start_a]),&(x[start_y]),len_y,&(y[start_x]),len_x,store,ierr);
                    }
                  pntr++;
                }
            }
        }
      *ierr = 0;
    }else  if((type=='H')&&(!(part=='B'))&&(m==n)) {
      if (part=='U') {
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  start_y=mat->BP1[i]+ofs;
                  len_y=mat->BP1[i+1]-mat->BP1[i];

                  //IA1 block colume id to block id
                  start_x = mat->BP2[mat->IA1[pntr+ofs]+ofs] + ofs;
                  len_x = mat->BP2[mat->IA1[pntr+ofs]+ofs+1]-mat->BP2[mat->IA1[pntr+ofs]+ofs];

                  start_a = mat->IA2[pntr+ofs] + ofs;


                  if(i==mat->IA1[pntr+ofs] + ofs)
                    dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                  else  if (i<mat->IA1[pntr+ofs] + ofs)
                    {
                      dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                      dblock_H_mult_vec (&(mat->A[start_a]),&(x[start_y]),len_y,&(y[start_x]),len_x,store,ierr);
                    }
                  pntr++;
                }
            }
        }else{
          for(i=0;i<mb;i++)
            {
              pntr = mat->PB[i];
              while(pntr<mat->PE[i])
                {
                  start_y=mat->BP1[i]+ofs;
                  len_y=mat->BP1[i+1]-mat->BP1[i];

                  //IA1 block colume id to block id
                  start_x = mat->BP2[mat->IA1[pntr+ofs]+ofs] + ofs;
                  len_x = mat->BP2[mat->IA1[pntr+ofs]+ofs+1]-mat->BP2[mat->IA1[pntr+ofs]+ofs];

                  start_a = mat->IA2[pntr+ofs] + ofs;


                  if(i==mat->IA1[pntr+ofs] + ofs)
                    dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                  else  if (i>mat->IA1[pntr+ofs] + ofs)
                    {
                      dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
                      dblock_H_mult_vec (&(mat->A[start_a]),&(x[start_y]),len_y,&(y[start_x]),len_x,store,ierr);
                    }
                  pntr++;
                }
            }
        }
      *ierr = 0;
    }
  else
    { //no symmetry
      for(i=0;i<mb;i++)
        {
          pntr = mat->PB[i];
          while(pntr<mat->PE[i])
            {
              start_y=mat->BP1[i]+ofs;
              len_y=mat->BP1[i+1]-mat->BP1[i];

              //IA1 block colume id to block id
              start_x = mat->BP2[mat->IA1[pntr+ofs]+ofs] + ofs;
              len_x = mat->BP2[mat->IA1[pntr+ofs]+ofs+1]-mat->BP2[mat->IA1[pntr+ofs]+ofs];

              start_a = mat->IA2[pntr+ofs] + ofs;
              dblock_mult_vec (&(mat->A[start_a]),&(x[start_x]),len_x,&(y[start_y]),len_y,store,ierr);
              pntr++;
            }
        }

      *ierr = 0;
    }
}// drmbv_vbr

