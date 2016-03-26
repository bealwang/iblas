#ifndef HASH_H
#define HASH_H
#include "blas_sparse_namedconstants.h"
#include <malloc.h>
// #include <iostream>
// **********************************************************************
//     Author : C. Voemel
//
//     Date of last modification : 3.1.2002
//
//     Description : A hash table for 'COO' and 'BCO' triangular solver
// **********************************************************************

typedef struct
{
  int jndx,val_pos;
  capsule* pntr;
}capsule;

typedef struct
{
  capsule* pntr;
}cappntr;

void setup_hash(int n,int* ierr);
void new_capsule_main(int indx,int jndx,int pos,int* ierr);
void print_hash();
void remove_hash(int* ierr);
void del_capsule(int nmb,int* ierr);




#endif // HASH_H
