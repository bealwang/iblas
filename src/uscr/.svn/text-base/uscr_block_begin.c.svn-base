#include "uscr/USCR_BEGIN.h"
d_matrix* duscr_block_begin (int Mb,int Nb,int k,int l)
{
  int m;
  d_matrix* dpmatrix;
  m=1;
  if((Mb<=0)||(Nb<=0)){// then
      return NULL;
    }
  else
    {
      dpmatrix=new_d_matrix(m);
      dpmatrix=daccess_matrix(dpmatrix);
      if(dpmatrix==NULL) return NULL;
      dpmatrix->DIM[3]=Mb;     //nb_of_block_rows
      dpmatrix->DIM[4]=Nb;     //nb_of_block_cols
      dpmatrix->DIM[5]=k;      //nb_of_rows_in_block
      dpmatrix->DIM[6]=l;      //nb_of_cols_in_block
      dpmatrix->format='b';
    }
  return dpmatrix;
}
