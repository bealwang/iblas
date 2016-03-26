#include "properties.h"


void get_descra(char* descra,int descriptor,char* message,int *ierr)
{
  *ierr =-1;
  //*message='';
  switch (descriptor)
    {
    case 'a':
      *message=descra[2];
      break;
    case 'b':
      *message=descra[4];
      break;
    case 'd':
      *message=descra[0];
      break;
    case 'f':
      *message=descra[6];
      break;
    case 'r':
      *message=descra[1];
      break;
    case's':
      *message=descra[3];
      break;
    case't':
      *message=descra[5];
      break;
    default:
      return;
    }
  *ierr =0;
}

void set_descra(char* descra, int prpty,int *ierr)
{
  int dummy;
  //*descra='';
  *ierr =-1;
  dummy=prpty;
  //check if matrix has unstroed unity diagomal
  if (dummy%2==1)
    descra[0]='U';
  else
    descra[0]='N';
  dummy=dummy - dummy%2;
  // repeated indices
  if (dummy%4 ==2)
    descra[1]='R';
  else
    descra[1]='U'; //default
  dummy =dummy - dummy % 4;
  //both/lower/upper half of matrix specified
  switch (dummy%16)
    {
    case  0:
      descra[2] ='B'; //default
      break;
    case 4:
      descra[2] ='U' ;
      break;
    case 8:
      descra[2] ='L' ;
      break;
    default:
      return ;
    }
      dummy =dummy - dummy %16;
      //matrix is irregular/reguarl/unassebmled
      switch (dummy % 64)
        {
        case 0:
          descra[3] ='I'; //default
          break;
        case 16:
          descra[3] ='R';
          break;
        case 32:
          descra[3] ='U';
          break;
        default:
          return ;
        }
      dummy =dummy - dummy % 64;
      //index base
      if (dummy%128==64)
        descra[4]='C';
      else
          descra[4]='F';

      dummy =dummy - dummy % 128;
      //matrix type
      switch (dummy % 1024)
        {
        case 0:
          descra[5] ='G';
          break;
        case 128:
          descra[5] ='S'; //default
          break;
        case 256:
          descra[5] ='H';
          break;
        case (512):
          descra[5] ='T';
          break;
        default:
          return;
        }
      dummy =dummy - dummy % 1024;
      //row major or column major
      if (dummy % 2048 ==1024)
        descra[6]='R';
      else
        descra[6]='C';

  *ierr =0;
}


void get_infoa(int* infoa,char descr, int* val,int *ierr)
{
  *val=-1;
  *ierr =-1;
  switch (descr)
    {
    case 'b': //base of array indices
      *val =infoa[1];
      break;
    case 'c':  //copy or not
      *val =infoa[8];
      break;
    case 'd':  //multidim array:row-dim of block
      *val =infoa[2];
      break;
    case 'e': //multidim arrau: col-dim of block
      *val =infoa[3];
      break;
    case 'f': //block structure array: row-dim inblocks
      *val =infoa[4];
      break;
    case 'g': //block structure array: col-dim inblocks
      *val =infoa[5];
      break;
    case 'n':
      *val =infoa[0];
      break;
    default:
      return ;
    }
  *ierr=0;
  return;
}

void set_infoa(int* infoa, char descr,int val,int *ierr)
{
  *ierr=-1;
  if(val<0) return;
  switch (descr)
    {
    case 'b': //base of array indices
      infoa[1] =val;
      break;
    case('c'): //copy or not
      infoa[8] =val;
      break;
    case('d'): //multidim array:row-dim of blocks
      infoa[2] =val;
      break;
    case('e'): //multidim array:col-dim of blocks
      infoa[3]=val;
      break;
    case('f'): //Block structure array:row-dim in blocks
      infoa[4] =val;
      break;
    case('g'): //Block structure array:col-dim in blocks
      infoa[5] =val;
      break;
    case('n'): //nnz
      infoa[0] =val;
      break;
    default:
      return;
    }
  *ierr=0;
}

