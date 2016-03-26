#include "gpu/gpu_rmbv_csc.h"
#include <string.h>
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS MV MULT. WITH MATRIX IN 'CSC'-STORAGE
//                   rmbv = Right Multiplication By Vector: y=Ax
void drmbv_csc (DSPMAT* mat, double* x, int n, double* y, int m, int* ierr)
{
  int base,ofs,j,pntr;
  char diag,type,part;
  *ierr=-1;
  if ((mat->FIDA!=CSC_FORMAT)||(mat->M!=m)||(mat->K!=n)) {
      *ierr=blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'b',&base,ierr);
 
  ofs=1 - base;
  ofs=0;// didn't konw its meaning defalut to 0
  get_descra(mat->DESCRA,'d',&diag,ierr);
  
  get_descra(mat->DESCRA,'t',&type,ierr);
  
  get_descra(mat->DESCRA,'a',&part,ierr);
  if (*ierr!=0) {
      *ierr=blas_error_param;
      return;
    }
  for(j=0;j<m;j++)
  y[j]=0;

  if (diag=='U') { //process unstored diagonal
      if (m==n) {//若有对角线，一般来说行列是相等的
          for(j=0;j<m;j++)
            y[j]=x[j];//because diagonal defalut to be 1
        }
      else
        {
          *ierr=blas_error_param;
          return;
        }
    }

  //type equal to symmetic
  if ((type=='S'||type=='H')&&(!(part=='B'))&&(m==n)) {//共轭矩阵只是针对复数而言,故对于double来说两者一样
      if (part=='U')//upper
        for(j=0;j<mat->K;j++)//j is index of column
          {
            pntr=mat->PB[j];
            while(pntr<mat->PE[j])
              {
                if(j==mat->IA1[pntr+ofs]+ofs)//行号等于列号的情况下，即对角线上
                  y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                else if (j>mat->IA1[pntr+ofs]+ofs)//上三角部分,因为是对称，所以加两边
                  {//右乘中，左矩阵值的列号等于右向量的行号，左矩阵的行号为结果向量的行号
                    y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                    y[j]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                  }
                pntr++;
              }
          }
      else//lower
        for(j=0;j<mat->K;j++)
          {
            pntr=mat->PB[j];
            while(pntr<mat->PE[j])
              {
                if(j==mat->IA1[pntr+ofs]+ofs)
                    y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                else if (j<mat->IA1[pntr+ofs]+ofs)
                  {
                    y[mat->IA1[pntr+ofs]+ofs]+=mat->A[pntr+ofs]*x[j];
                    y[j]+=mat->A[pntr+ofs]*x[mat->IA1[pntr+ofs]+ofs];
                  }
                pntr++;
              }
          }
      *ierr=0;
    }
      else{// unsymmtric matrix, Yi=Aij*Xj,i=0,...k
          #pragma omp parallel for num_threads(dtn(mat->K,MIN_ITERATOR_NUM))
          for(j=0;j<mat->K;j++)
            {
              int col_begin = mat->PB[j];
              int col_end = mat->PE[j];
              int i;
              for(i=col_begin;i<col_end;i++)
                {
                  y[mat->IA1[i]]+=mat->A[i]*x[j];
                }
            }

          // rmbv_csc(mat->IA1,mat->A,mat->PB,mat->PE,mat->M,mat->K,mat->n_IA1,x,y);
          *ierr=0;
        }

}
