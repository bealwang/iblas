#include "uscr/usds.h"
int BLAS_usds(void* sp_l,int type)
{
int ierr;

DSPMAT* dsp_data;
int val;
switch(type){
  case 3://DSP_MATRIX:
    // **********************************************************************
    dsp_data=accessdata_dsp(sp_l,&ierr);
    if (dsp_data==NULL){// then
        ierr = blas_error_param;
        return ierr;
      }//end if
    get_infoa(dsp_data->INFOA,'c',&val,&ierr);
    if (ierr!=0) {//then
        ierr = blas_error_param;
        return ierr;
      }//end if
    if(val==COP_OF_SOURCE)
      {//then
        // *** Deallocate extra storage for copy of matrix *** //

        if(dsp_data->FIDA==COO_FORMAT||dsp_data->FIDA==BCO_FORMAT){
            aligned_free(dsp_data->A);
            aligned_free(dsp_data->IA1);
            aligned_free(dsp_data->IA2);
          }
        else if(dsp_data->FIDA==CSC_FORMAT||dsp_data->FIDA==BSC_FORMAT){

            aligned_free(dsp_data->A);
            aligned_free(dsp_data->IA1);
            aligned_free(dsp_data->PB);
            aligned_free(dsp_data->PE);
          }
        else if(dsp_data->FIDA==CSR_FORMAT||dsp_data->FIDA==BSR_FORMAT){
            aligned_free(dsp_data->A);
            aligned_free(dsp_data->IA1);
            aligned_free(dsp_data->PB);
            aligned_free(dsp_data->PE);
          }
        else if(dsp_data->FIDA==DIA_FORMAT ||dsp_data->FIDA==BDI_FORMAT){
            aligned_free(dsp_data->A);
            aligned_free(dsp_data->IA1);
          }
        else if(dsp_data->FIDA==VBR_FORMAT)
          {
            aligned_free(dsp_data->A);
            aligned_free(dsp_data->IA1);
            aligned_free(dsp_data->IA2);
            aligned_free(dsp_data->PB);
            aligned_free(dsp_data->PE);
            aligned_free(dsp_data->BP1);
            aligned_free(dsp_data->BP2);
          }
        else
          {
            ierr = blas_error_param;
            return ierr;
          }
        if(ierr!=0) {//then
            ierr=blas_error_memdeloc;
            return ierr;
          }//end if
      }//end if
    del_dsp (sp_l,&ierr);
    if (ierr!=0){//then
        ierr=blas_error_memdeloc;
        return ierr;
      }

    break;


}
if (ierr!=0){//then
  ierr = blas_error_param;
  return ierr ;
}
return ierr;
//end if
}//end subroutine usds
