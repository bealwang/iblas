#include <assert.h>
#include "INSERTING.h"

d_matrix* new_d_matrix (int Mb)
{
  d_matrix* d_matrix_start;
  int i=0;
 //if (d_matrix_start==NULL) {//then
      d_matrix_start=(d_matrix*)aligned_malloc(sizeof(d_matrix));
      d_matrix_start->number=DSP_MATRIX;
      d_matrix_start->number=-d_matrix_start->number;
      d_matrix_start->pntr=NULL;
   /* }else{//else
      matrix_insert=malloc(sizeof(d_matrix));
      matrix_insert->number=d_matrix_start->number-no_of_types;
      matrix_insert->pntr=d_matrix_start ;
      d_matrix_start=matrix_insert;
    }*///end if

      for(i=0;i<7;i++)
        d_matrix_start->DIM[i]=0;

  d_matrix_start->property=blas_general+blas_zero_base+blas_col_major;
  d_matrix_start->news=1 ;   //new=0:blas_open_handle, new=1: blas_new_handle
  d_matrix_start->format=' ';
  d_matrix_start->sub_cols=NULL;
  d_matrix_start->sub_rows=NULL;
  d_matrix_start->d_element_start=NULL;
  d_matrix_start->trb=(int*)aligned_malloc(sizeof(int)*Mb);
  d_matrix_start->tre=(int*)aligned_malloc(sizeof(int)*Mb);
  return d_matrix_start;
}//end subroutine new_i_matrix

void dealloc_d_matrix  (d_matrix* matrix)
{
  assert((matrix->trb!=NULL) && (matrix->tre!=NULL));
  aligned_free(matrix->trb);
  aligned_free(matrix->tre);
  aligned_free(matrix);
}//end subroutine dealloc_d_matrix

d_matrix* daccess_matrix (d_matrix* p)
{
  d_matrix* pmatrix;
  d_matrix* matrix_tester;
  return p;

}//end subroutine iaccess_matrix

int new_d_element (d_matrix* pmatrix,int* istat)
{
  int nmb_element;
  d_element* element_insert;
  *istat=-1;
  if (pmatrix->d_element_start==NULL) {//then
      pmatrix->news=0 ;//status changed to blas_open_handle
      pmatrix->d_element_start=(d_element*)aligned_malloc (sizeof(d_element));
      pmatrix->d_element_start->number=1; //will certainly changed
      pmatrix->d_element_start->pntr=NULL;
    }else{//else
      element_insert=(d_element*)aligned_malloc(sizeof(d_element));
      element_insert->pntr=pmatrix->d_element_start;
      element_insert->number=pmatrix->d_element_start->number+1;
      pmatrix->d_element_start=element_insert;
    }//end if

  switch (pmatrix->format)
    {
    case 'n':
      pmatrix->d_element_start->contents.pntin.value=0;
      pmatrix->d_element_start->contents.pntin.row_ind=-1;
      pmatrix->d_element_start->contents.pntin.col_ind=-1;
      pmatrix->d_element_start->contents.blin.value=NULL;
      pmatrix->d_element_start->contents.vblin.value=NULL;
      break;
    case 'b':
      pmatrix->d_element_start->contents.blin.value=NULL;
      pmatrix->d_element_start->contents.vblin.value=NULL;
      pmatrix->d_element_start->contents.blin.row_block_ind=-1;
      pmatrix->d_element_start->contents.blin.col_block_ind=-1;
      break;
    case 'v':
      pmatrix->d_element_start->contents.blin.value=NULL;
      pmatrix->d_element_start->contents.vblin.value=NULL;
      pmatrix->d_element_start->contents.vblin.row_vblock_ind=-1;
      pmatrix->d_element_start->contents.vblin.col_vblock_ind=-1;
      break;
    default:
      *istat=blas_error_param;
      return -1;
    }
  nmb_element=pmatrix->d_element_start->number;
  *istat=0;
  return nmb_element;
}//end subroutine new_d_element

void dealloc_d_element  (int nmb_element,d_matrix* pmatrix,int* istat)
{
  d_element* element_tester;
  *istat=-1;
  if(pmatrix->d_element_start->pntr==NULL) {//then
      if(pmatrix->d_element_start->number==nmb_element) {//then
          if(pmatrix->d_element_start->contents.vblin.value!=NULL)
            aligned_free(pmatrix->d_element_start->contents.vblin.value);

          if(pmatrix->d_element_start->contents.blin.value!=NULL)
            aligned_free(pmatrix->d_element_start->contents.blin.value);

          if (pmatrix->d_element_start!=NULL)
            aligned_free(pmatrix->d_element_start);
          pmatrix->d_element_start==NULL;
        }//end if
      *istat=0;
      return;
    }
  else
    {
      element_tester=pmatrix->d_element_start;
      if(element_tester->number==nmb_element) {//then
          pmatrix->d_element_start=element_tester->pntr;

          if(element_tester->contents.vblin.value!=NULL)
              aligned_free(element_tester->contents.vblin.value);

          if(element_tester->contents.blin.value!=NULL)
              aligned_free(element_tester->contents.blin.value);

          if (element_tester!=NULL)
              aligned_free(element_tester);
          *istat=0;
          return;
        }//end if
      element_tester=pmatrix->d_element_start->pntr;
      while(element_tester!=NULL){//while
          if(element_tester->number==nmb_element){//then
              if(element_tester->contents.vblin.value!=NULL)
                aligned_free(element_tester->contents.vblin.value);

              if(element_tester->contents.blin.value!=NULL)
                aligned_free(element_tester->contents.blin.value);

              if (element_tester!=NULL)
                aligned_free(element_tester);

              istat=0;
              return;
            }
          else
            {
              element_tester=element_tester->pntr;
            }//end if
        }//end do
    }//end if
  }//end subroutine dealloc_d_element

d_inelement* daccess_element (int nmb_element,d_matrix* pmatrix)
{
  d_inelement* pelement;
  d_element* element_tester;
  pelement=NULL;
  element_tester=pmatrix->d_element_start;
  while((element_tester->number!=nmb_element)&&(element_tester->pntr!=NULL))
    element_tester=element_tester->pntr;

  if (element_tester->number==nmb_element) {//then
      pelement=&(element_tester->contents);
    }else{//else
      pelement=NULL;
    }//end if
  return pelement;
}//end subroutine  daccess_element

//find the element(i,j) and return its nmb
int d_element_num  (d_matrix* pmatrix,int i,int j,int* istat)
{

  int nmb_element;
  d_element* element_tester;
  int finder;
  *istat=-1;
  switch (pmatrix->format){
    case 'n':
      element_tester=pmatrix->d_element_start;
      if(element_tester->pntr==NULL) {//then
          if((element_tester->contents.pntin.row_ind==i)&&(element_tester->contents.pntin.col_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }else{
          finder=FALSE;
          while((element_tester->pntr!=NULL)&&(!finder)){//while
              if((element_tester->contents.pntin.row_ind==i)&&(element_tester->contents.pntin.col_ind==j)){//then
                  finder=TRUE;
                }else{
                  element_tester=element_tester->pntr;
                }//end if
            }//end do
          if((element_tester->contents.pntin.row_ind==i)&&(element_tester->contents.pntin.col_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }//end if
      break;
    case 'b':
      element_tester=pmatrix->d_element_start;
      if(element_tester->pntr==NULL) {//then
          if((element_tester->contents.blin.row_block_ind==i)&&(element_tester->contents.blin.col_block_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }else{
          finder=FALSE;
          while((element_tester->pntr!=NULL)&&(!finder)){//while
              if((element_tester->contents.blin.row_block_ind==i)&&(element_tester->contents.blin.col_block_ind==j)){//then
                  finder=TRUE;
                }else{
                  element_tester=element_tester->pntr;
                }// end if
            }//end do
          if((element_tester->contents.blin.row_block_ind==i)&&(element_tester->contents.blin.col_block_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }//end if
      break;
    case 'v':
      element_tester=pmatrix->d_element_start;
      if(element_tester->pntr==NULL){//then
          if((element_tester->contents.vblin.row_vblock_ind==i)&&(element_tester->contents.vblin.col_vblock_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }else{
          finder=FALSE;
          while((element_tester->pntr!=NULL)&&(!finder)){
              if((element_tester->contents.vblin.row_vblock_ind==i)&&(element_tester->contents.vblin.col_vblock_ind==j)){//then
                  finder=TRUE;
                }else{
                  element_tester=element_tester->pntr;
                }//end if
            }//end do
          if((element_tester->contents.vblin.row_vblock_ind==i)&&(element_tester->contents.vblin.col_vblock_ind==j)){//then
              nmb_element=element_tester->number;
            }else{
              nmb_element=0;
            }//end if
        }//end if
      break;
    default:
      *istat=blas_error_param;
      return -1;
    }
  *istat=0;
  return nmb_element;
}//end subroutine d_element_num

void d_dealloc (d_matrix* pmatrix,int* istat)
{
  d_element *element_tester,*next_element;
  *istat=-1;
  pmatrix=daccess_matrix(pmatrix);
  element_tester=pmatrix->d_element_start;
  if(element_tester->pntr==NULL){//then
      dealloc_d_element(element_tester->number,pmatrix,istat);
      if (*istat!=0) return;
    }
  else
    {
      next_element=element_tester->pntr;
      while(next_element!=NULL){//while
          dealloc_d_element (element_tester->number,pmatrix,istat);
          if (*istat!=0) return;
          element_tester=next_element;
          next_element=element_tester->pntr;
        }
      dealloc_d_element(element_tester->number,pmatrix,istat);
 if (*istat!=0) return;
    }
  dealloc_d_matrix (pmatrix);

  *istat=0;
  return;
}//end subroutine d_dealloc
