#include <stdio.h>
#include "hash.h"
// **********************************************************************
//     Author : C. Voemel
//
//     Date of last modification : 3.1.2002
//
//     Description : A hash table for 'COO' and 'BCO' triangular solver
// **********************************************************************
capsule* hash;
cappntr* hash_top;
int hash_length;

void setup_hash(int n,int* ierr){
  int i;
  *ierr=-1;
  hash=NULL;
  hash=(capsule* )malloc (sizeof(capsule)*n);
  hash_top=(cappntr *)malloc (sizeof(cappntr)*n);
  hash_length=n;

  if (hash_top==NULL||hash==NULL)
    {
      *ierr=blas_error_memalloc;
      return;
    }
  for(i=0;i<n;i++)
    {
      hash[i]->pntr=NULL;
      hash->jndx= -1;
      hash->val_pos= -1;
      hash_top[i]->pntr=hash[i];
    }
  *ierr=0;
}

void new_capsule_main(int indx,int jndx,int pos,int* ierr)
{
  capsule* cap=NULL;
  *ierr = -1;
  //判断indx是否越界
  if ((indx<0)||(indx>=hash_length))
    return;

  if(indx==jndx)
    {
      hash[indx]->val_pos = pos;
      hash[indx]->jndx = jndx;
    }
  else
    {
      cap=(capsule *)malloc(sizeof(capsule));
      if (cap==NULL)
        {
          *ierr=blas_error_memalloc;
          return;
        }
      cap->val_pos = pos;
      cap->jndx = jndx;
      cap->pntr=NULL;
      hash_top[indx]->pntr->pntr =cap;
      hash_top[indx]->pntr =cap;
    }
  *ierr = 0;
}


void print_hash(){
  int i;
  capsule* dummy=NULL;
  for( i=0;i<hash_length;i++)
    {
      printf("print_hash(%d) ",i);
      dummy=hash[i];
      while(dummy->pntr!=NULL)
        {
          printf("jndx : %d",dummy->jndx);
          printf("val_pos : %d",dummy->val_pos);
          dummy=dummy->pntr;
        }
      printf("jndx :%d",dummy->jndx);
      printf("val_pos : %d",dummy->val_pos);
    }
}

void remove_hash(int* ierr)
{
  int i;
  *ierr = -1;
  for( i=0;i<hash_length;i++)
    while(hash_top[i]->pntr!=hash[i])
      {
       del_capsule(i,ierr);
        if (*ierr!=0)
          {
            *ierr=blas_error_memdeloc;
            return;
          }
      }

  delete hash;
  delete hash_top;
  hash=NULL;
  hash_top=NULL;

  if (*ierr!=0)
    {
      *ierr=blas_error_memdeloc;;
      return;
    }
}



void del_capsule(int nmb,int* ierr)
{
  capsule* dummy=NULL;
  dummy=hash[nmb];
  if (dummy==hash_top[nmb]->pntr)
    {
      *ierr = -1;
      return;
    }

  while((dummy->pntr)!=(hash_top[nmb]->pntr))
    dummy =dummy->pntr;

  hash_top[nmb]->pntr=dummy;

  delete dummy->pntr;
  dummy->pntr=NULL;

  if(*ierr!=0)
    {
      *ierr=blas_error_memdeloc;
      return;
    }

}
