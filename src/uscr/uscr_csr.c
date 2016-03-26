#include "uscr/uscr_csr.h"
dsp_linknode* duscr_csr (int m,int n,double* val,int n_val,int* indx,int n_indx,int* pntrb,int n_pntrb,int* pntre,int n_pntre,int prpty,int* istat)
{
  int options,base,nnz;
  bool  COPY;
  char  message;
  DSPMAT* dsp_data;
  dsp_linknode* dsp_l;
  int i;
  options=*istat;
  *istat=-1; //if not changed later,routine has failed
  dsp_l=new_dsp(istat);
  if(*istat!=0) {//then
      *istat=blas_error_memalloc;
      return NULL;
    }//end if
  dsp_data=accessdata_dsp(dsp_l,istat);
  if(*istat!=0) {//then
      *istat=blas_error_param;
      return NULL;
    }//end if
  dsp_data->FIDA=CSR_FORMAT;
  dsp_data->M=m;
  dsp_data->K=n;
  set_descra(dsp_data->DESCRA,prpty,istat);
  get_descra(dsp_data->DESCRA,'b',&message,istat);

  if (message=='C')
      base=C_BASE;
  else //Assuming F base
    base=F_BASE;

  set_infoa(dsp_data->INFOA,'b',base,istat);
  nnz=maxval(pntre,n_pntre)-base;
  set_infoa(dsp_data->INFOA,'n',nnz,istat);
  if((nnz!=n_indx)||(m!=n_pntrb)||(minval(indx,n_indx)<base)||(maxval(indx,n_indx)>n-1+base)||(m!=n_pntre)||(nnz!=n_val)) {//then
      BLAS_usds(dsp_l,3);
      *istat=blas_error_param;
      return NULL;
    }

  //init the size of array in dsp
  dsp_data->n_A=n_val;
  dsp_data->n_IA1=n_indx;
  dsp_data->n_PB=n_pntrb;
  dsp_data->n_PE=n_pntre;
  dsp_data->n_BP1=0;
  dsp_data->n_BP2=0;
  dsp_data->n_IA2=0;


  if(options>0) {//then
      // decision rule whether or not to copy
      COPY=TRUE;
      if(COPY) {//then
          options=-1; //copy
        }else{
          options=0;  //reference
        }//end if
    }//end if
  if (options==0) {//then
      set_infoa(dsp_data->INFOA,'c',REF_OF_SOURCE,istat);
      if (*istat!=0) {
          *istat=blas_error_param;
          return NULL;
        }

      // create reference to original matrix
      dsp_data->A=val;
      dsp_data->IA1=indx;
      dsp_data->PB=pntrb;
      dsp_data->PE=pntre;

      *istat=0;
    }else{
      // The additional required memory is DEALLOCATED later in USDS//
      set_infoa(dsp_data->INFOA,'c',COP_OF_SOURCE,istat);
      if (*istat!=0) {//
          *istat=blas_error_param;
          return NULL;
        }
      // copy original data
      dsp_data->A=(double*)aligned_malloc(sizeof(double)*n_val);
      dsp_data->IA1=(int*)aligned_malloc(sizeof(int)*n_indx);
      dsp_data->PB=(int*)aligned_malloc(sizeof(int)*n_pntrb);
      dsp_data->PE=(int*)aligned_malloc(sizeof(int)*n_pntre);

      if (*istat!=0) {//
          *istat=blas_error_memalloc;
          return NULL;
        }
      for(i=0;i<dsp_data->n_A;i++)
        dsp_data->A[i]=val[i];
      for(i=0;i<dsp_data->n_IA1;i++)
        dsp_data->IA1[i]=indx[i];
      for(i=0;i<dsp_data->n_PB;i++)
        dsp_data->PB[i]=pntrb[i];
      for(i=0;i<dsp_data->n_PE;i++)
        dsp_data->PE[i]=pntre[i];
      *istat=1;
    }


  if(*istat>=0)
    {
      *istat=0;
      return dsp_l;
    }
  else
    return NULL;
}
