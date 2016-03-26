#ifndef INSERTING_H
#define INSERTING_H
#include "properties.h"
#include "comm_tools.h"
#include <malloc.h>
typedef struct
{
  int row_ind,col_ind;
  int value;
}i_inpnt1;

typedef struct{
  int row_block_ind,col_block_ind;
  int** value;
}i_inblock;

typedef struct {
  int row_vblock_ind,col_vblock_ind;
  int** value;
}i_invblock;

typedef struct{
  i_inblock blin;
  i_inpnt1 pntin;
  i_invblock vblin;
}i_inelement;

typedef struct ielement{
  int number;
  i_inelement contents;
  struct ielement* pntr;
}i_element;

typedef struct imatrix{
  int DIM[7];
  int property,number,news;
  char format;//n,b,v means:normal,block,vblock
  int *sub_rows,*sub_cols,*trb,*tre;
  i_element* i_element_start;
  struct imatrix* pntr;
}i_matrix;

typedef struct
{
  int row_ind,col_ind;
  double value;
}d_inpnt1;

typedef struct{
  int row_block_ind,col_block_ind;
  double** value;
}d_inblock;

typedef struct {
  int row_vblock_ind,col_vblock_ind;
  double** value;
}d_invblock;

typedef struct{
  d_inblock blin;
  d_inpnt1 pntin;
  d_invblock vblin;
}d_inelement;

typedef struct delement{
  int number;
  d_inelement contents;
  struct delement* pntr;
}d_element;

typedef struct dmatrix{
  int DIM[7];
  int property,number,news;
  char format;
  int *sub_rows,*sub_cols,*trb,*tre;
  d_element* d_element_start;
  struct dmatrix* pntr;
}d_matrix;

d_inelement* daccess_element (int nmb_element,d_matrix* pmatrix);

d_matrix* new_d_matrix (int Mb);

void dealloc_d_matrix  (d_matrix* matrix);

d_matrix* daccess_matrix (d_matrix* p);

int new_d_element (d_matrix* pmatrix,int* istat);

void dealloc_d_element  (int nmb_element,d_matrix* pmatrix,int* istat);

d_inelement* daccess_element (int nmb_element,d_matrix* pmatrix);

//find the element(i,j) and return its nmb
int d_element_num  (d_matrix* pmatrix,int i,int j,int* istat);

void d_dealloc (d_matrix* pmatrix,int* istat);

#endif
