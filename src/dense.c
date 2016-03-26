#include "dense.h"
#include "blas_sparse_namedconstants.h"
/*
*A is a two dimmenision array stored in one dimension
*A is n*m
*/
#include "dense.h"
//done
void dblock_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr)
{
  *ierr=-1;
  if (store=='C')//column major use right mutil
      dblock_r_mult_vec(A,x,n,y,m,ierr);
  else
      dblock_l_mult_vec(A,x,n,y,m,ierr);
}
//done
void dblock_Z_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr)
{
  *ierr=-1;
  if (store=='C')//column major use right mutil
      dblock_r_mult_vec(A,x,n,y,m,ierr);
  else
      dblock_l_mult_vec(A,x,n,y,m,ierr);
   //y=(y)
}

//transform 对称变化乘
void dblock_T_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr)
{
  *ierr=-1;
  if (store=='C')//column major use right mutil
      dblock_l_mult_vec(A,x,n,y,m,ierr);
  else
      dblock_r_mult_vec(A,x,n,y,m,ierr);
}

//H 共轭变化乘  是Z和T的混合，但对double型无效
void dblock_H_mult_vec (double *A, double *x, int n, double *y, int m, char store, int* ierr)
{
  *ierr=-1;
   if (store=='C')//column major use right mutil
       dblock_l_mult_vec(A,x,n,y,m,ierr);
   else
       dblock_r_mult_vec(A,x,n,y,m,ierr);
   //y=(y)
}

//done:y=xA
void dblock_r_mult_vec (double *A, double *x, int n, double *y, int m,int* ierr)
{
  int j;
  *ierr=-1;
  printf("R\n");
  for (j=0;j<n;j++){
    int i;
    for (i=0;i<m;i++)
      {
        y[i]+=A[j*m+i]*x[j];
      }
  }
  *ierr=0;
}

//done:y=Ax
void dblock_l_mult_vec (double *A, double *x, int n, double *y, int m,int* ierr)
{   int j,i;
    *ierr=-1;
    printf("L\n");
    for (j=0;j<m;j++)
       for(i=0;i<n;i++)
          y[j]+=A[j*n+i]*x[i];
    *ierr=0;
}










/*
//dinvert
void dinvert_left_lower (double *A, double *x, int n, double *y, int m, char store, int* ierr)
{  if (store=='C')
    {
      dinvert_r_left_lower(A,x,n,y,m,ierr);
    }
  else
    {
      dinvert_l_right_upper(A,x,n,y,m,ierr);
    }
   return 0;
}

//
void dinvert_T_left_lower (double *A, double *x, int n, double *y, int m, char store,
                      int* ierr)
{  if (store=='C')
    {
      dinvert_l_left_lower(A,x,n,y,m,ierr);
    }
  else
    {      dinvert_r_right_upper(A,x,n,y,m,ierr);
    }
   return 0;
}

//
void dinvert_T_right_upper (double *A, double *x, int n, double *y, int m, char store,int* ierr)
{  if (store=='C')
    {
      dinvert_l_right_upper(A,x,n,y,m,ierr);
    }
  else
    {      dinvert_r_left_lower(A,x,n,y,m,ierr);
    }
   return 0;
}

//left_lower,stored column-wise
//AY=x
//compute Y,and final output the result Y to input array x
//I think here left_lower means the left_lower is zero
void dinvert_r_left_lower (double *A, double *x, int n, double *y, int m, char store,int* ierr)
{  int j;
   for (j=0;j<n;j++)
     {
       if (A[j*n+j]!=0)	//make sure diagonal is not zero
         x[j]=x[j]/A[j*n+j];
       else
         {
           *ierr=blas_error_singtria;
           return;
         }
           for (i=j+1;i<n;i++)
           x[i]-=x[j]*A[j*n+i];

     }
    return 0;
}

//y=x*A^T
void dinvert_l_left_lower (double *A, double *x, int n, double *y, int m, char store,
                      int* ierr)
{  int j;
   for (j=0;j<n;j++)
     {      x[j]=x[j] - dot_product (A[j * n], x[j], j);
       if (A[(j - 1) * n + j] !=0)
         x[j]=x[j] / A[(j - 1) * n + j];
     }
    return 0;
}

void dinvert_r_right_upper (double *A, double *x, int n, double *y, int m, char store,
                       int* ierr)
{  int j;
   for (j=n;j >=0;j--)
     {      if (A[(j - 1) * n + j] !=0)
         x[j]=x[j] / A[(j - 1) * n + j];
       for (i=0;i<j;i++)
         {          x[i]=x[i] - x[j] * A[j * n + i];
         }
     }
    return 0;
}

void dinvert_l_right_upper (double *A, double *x, int n, double *y, int m, char store,
                       int* ierr)
{  int j;
   for (j=n - 1;j >=0;j--)
     {      x[j]=x[j] - dot_product (A[j * n], x[j], n - j);
       if (A[(j - 1) * n + j] !=0)
         x[j]=x[j] / A[(j - 1) * n + j];
     }
    return 0;
}
*/
