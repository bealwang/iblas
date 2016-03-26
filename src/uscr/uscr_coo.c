#include "uscr/uscr_coo.h"

dsp_linknode* duscr_coo(int m,int n,double* val,int n_val,int* indx,int* jndx,int nnz,int prpty,int* istat)
{
  int  ierr,options,base;
  bool  COPY;
  char  message;
  DSPMAT* dsp_data;
  dsp_linknode* dsp_l;
  int i;
  options=*istat;
  *istat=-1; //if not changed later,routine has failed
  ierr=-1;
  
  dsp_l=new_dsp(&ierr);
  if(ierr!=0) {//then
      ierr=blas_error_memalloc;
      return NULL;
    }//end if
  dsp_data=accessdata_dsp(dsp_l,&ierr);
  if(ierr!=0) {//then
      ierr=blas_error_param;
      return NULL;
    }//end if
  dsp_data->FIDA=COO_FORMAT;
  dsp_data->M=m;
  dsp_data->K=n;
  set_descra(dsp_data->DESCRA,prpty,&ierr);
  dsp_data->spc.M = m;
  dsp_data->spc.N = n;
  dsp_data->spc.NNZ = n_val;
  get_descra(dsp_data->DESCRA,'b',&message,&ierr);
  if(ierr!=0) {//then
      ierr=blas_error_param;
      return NULL;
    }//end if
  if(message=='C') {//then
      base=C_BASE;
    }else{ //Assuming F base
      base=F_BASE;
    }//end if
  set_infoa(dsp_data->INFOA,'b',base,&ierr);
  set_infoa(dsp_data->INFOA,'n',nnz,&ierr);
  if(ierr!=0) {//then
      ierr=blas_error_param;
      return NULL;
    }//end if
  if((nnz!=n_val)||(minval(indx,n_val)<base)||(minval(jndx,n_val)<base)||(maxval(indx,n_val)>m-1+base)||(maxval(jndx,n_val)>n-1+base)) {//then
      BLAS_usds(dsp_l,3);
      ierr=blas_error_param;
      return NULL;
    }//end if
  if(options>0) {//then
      // decision rule whether or not to copy
      COPY=TRUE;
      if(COPY) {//then
          options=-1; //copy
        }else{
          options=0;  //reference
        }//end if
    }//end if
  dsp_data->n_A=n_val;
  dsp_data->n_IA1=n_val;
  dsp_data->n_IA2=n_val;



  if(options==0) {//then
      set_infoa(dsp_data->INFOA,'c',REF_OF_SOURCE,&ierr);
      if(ierr!=0) {//then
          ierr=blas_error_param;
          return NULL;
        }//end if
      // create reference to original matrix
      dsp_data->A=val;
      dsp_data->IA1=indx;
      dsp_data->IA2=jndx;
      *istat=0;
    }else{
      // The additional required memory is DEALLOCATED later in USDS//
      set_infoa(dsp_data->INFOA,'c',COP_OF_SOURCE,&ierr);
      if(ierr!=0) {//then
          ierr=blas_error_param;
          return NULL;
        }//end if
      // copy original data
      dsp_data->A=(double*)aligned_malloc(sizeof(double)*n_val);
      dsp_data->IA1=(int*)aligned_malloc(sizeof(int)*n_val);
      dsp_data->IA2=(int*)aligned_malloc(sizeof(int)*n_val);

      //#pragma omp parallel for num_threads(dtn(n_val,MIN_ITERATOR_NUM))
      for(i=0;i<n_val;i++)
        {
          dsp_data->A[i]=val[i];
          dsp_data->IA1[i]=indx[i];
          dsp_data->IA2[i]=jndx[i];
        }
      *istat=1;
    }//end if
  if(*istat>=0)
    {
      *istat=0;
      return dsp_l;
    }
  else
    return NULL;
}//end subroutine duscr_coo
