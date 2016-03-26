#include "uscr/USCR_INSERT_CLIQUE.h"
#include "uscr/USCR_INSERT_ENTRY.h"
int BLAS_duscr_insert_clique (d_matrix* A, const int k, const int l, const double** val, const int row_stride, const int col_stride, const int *indx, const int *jndx)
{

  //      integer ::i,j,s_row,s_col
  int i,j;
  int istat=-1;
  
  //#pragma omp parallel for num_threads(dtn(l,MIN_ITERATOR_NUM)) private(i)
  for(j=0;j<l;j++)
    {  
      //int flag = 0;
      for(i=0;i<k;i++)
        {
          istat=BLAS_duscr_insert_entry(A,val[i][j],indx[i],jndx[j]);
          if(istat!=0){
            return istat;
            // flag = 1;
            // break;
          } 
        }
      // if(flag == 1){
      //   break;
      // }
    }
    return istat;
}
