// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 9.2.2012
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'COO'-STORAGE
//                   rsbv = Right Solve By Vector
// **********************************************************************
#include "hash.h"
#include "link.h"
#include  "properties.h"
//interface rsbv_coo
//  module procedure irsbv_coo
//  module procedure srsbv_coo
//  module procedure drsbv_coo
//  module procedure crsbv_coo
//  module procedure zrsbv_coo
//end interface
// **********************************************************************
// **********************************************************************
void irsbv_coo(ISPMAT* mat,int* x,int n_x,int &ierr)
{
  int i,n,base,ofs,nnz;
  char diag,part;
  capsule* dummy;
  ierr = -1;
  n = n_x;
  if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n))
    {
      ierr = blas_error_param;
      return;
    }
  get_infoa(mat->INFOA,'b',base,ierr);
  if (ierr!=0)
    {
      ierr = blas_error_param;
      return;
    }
  ofs = 1 - base;
  get_infoa(mat->INFOA,'n',nnz,ierr);
  if (ierr!=0)
    {
      ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'d',diag,ierr);
  if (ierr!=0)
    {
      ierr = blas_error_param;
      return;
    }
  get_descra(mat->DESCRA,'a',part,ierr);
  if (ierr!=0)
    {
      ierr = blas_error_param;
      return;
    }
  if ((part!='U')&&(part!='L'))
    {
      ierr = blas_error_param;
      return;
    }
  setup_hash(n,ierr);
  if (ierr!=0)
    {
      return;
    }
  for(i = 0;i<nnz;i++)
    {
      new_capsule_main(mat->IA1[i]+ofs,mat->IA2[i]+ofs,i,ierr);
      if (ierr!=0)
        {
          return;
        }
    }
  if (part=='L')
    {
      for(i=0;i<n;i++)
        {
          dummy = hash[i];
          while(dummy->pntr!=NULL)
            {
              dummy = dummy->pntr;
              x[i]=x[i]-x[dummy->jndx]*(mat->A[dummy->val_pos]);
            }
          if (diag!='U')
            {
              if(hash[i]->jndx==-1) {
                  ierr = blas_error_singtria;
                  return;
                }
              else
                {
                  if (mat->A[hash[i]->val_pos]!=0) {
                      x[i] = x[i]/(mat->A[hash[i]->val_pos]);
                    }
                  else
                    {
                      ierr = blas_error_singtria;
                      return;
                    }
                }
            }
        }
      ierr = 0;
    }
  else
    {
      for(i=n-1;i>=0;i--)
        {
          dummy = hash[i];
          while(dummy->pntr!=NULL)
            {
              dummy = dummy->pntr;
              x[i]=x[i]-x[dummy->jndx]*(mat->A[dummy->val_pos]);
            }
          if (diag!='U') {
              if(hash[i]->jndx==-1)
                {
                  ierr = blas_error_singtria;
                  return;
                }
              else
                {
                  if (mat->A[hash[i]->val_pos]!=0)
                    {
                      x[i] = x[i]/(mat->A[hash[i]->val_pos]);
                    }
                  else
                    {
                      ierr = blas_error_singtria;
                      return;
                    }
                }
            }
        }
      ierr = 0;
    }
  remove_hash(ierr);
}




