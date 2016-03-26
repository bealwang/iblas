#include "properties.h"
#include "conv/conv_tools.h"
#include "link.h"
#include "conv/usconv.h"
void usconv_dia2coo(void* spmat,int* ierr)
{
  int res,rest,NNZ,LDA,i;
  ISPMAT*  isp_data;
  SSPMAT*  ssp_data;
  DSPMAT*  dsp_data;
  CSPMAT*  csp_data;
  ZSPMAT*  zsp_data;
  (*ierr)=-1;
  rest=((dsp_linknode*)spmat)->number;
  switch (rest){
    case ISP_MATRIX:
      // **********************************************************************;

      break;
      // **********************************************************************;
    case SSP_MATRIX:
      // **********************************************************************;

      break;
      // **********************************************************************;
    case DSP_MATRIX:
      // **********************************************************************;
      dsp_data=accessdata_dsp(spmat,ierr);
      if (*ierr!=0) {
          *ierr = blas_error_param;
          return;
        }
      if(dsp_data->FIDA==DIA_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dsp_data->FIDA=COO_FORMAT;
              get_infoa(dsp_data->INFOA ,'d',&LDA,ierr);// !row-dim of val
              get_infoa(dsp_data->INFOA ,'n',&NNZ,ierr);//
              dpre_usconv_dia2coo (&(dsp_data->A), &(dsp_data->IA1),dsp_data->n_IA1 ,&(dsp_data->IA2),LDA,NNZ,dsp_data->M);
              dsp_data->n_A=NNZ;
              dsp_data->n_IA1=NNZ;
              dsp_data->n_IA2=NNZ;
              dsp_data->n_BP1=0;
              dsp_data->n_BP2=0;
              dsp_data->n_PB=0;
              dsp_data->n_PE=0;

              // printf("\n----------------DIA2COO:----------------\n");
              // for(i=0;i<dsp_data->n_IA1;i++){
              //   printf("IA1:%d  IA2:%d\n",dsp_data->IA1[i],dsp_data->IA2[i]);
              // }
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   printf("%lf  ",dsp_data->A[i]);
              // }
              // printf("\n----------------DIA2COO:----------------\n");

            }
        }else{
          *ierr = blas_error_param;
          return;
        }
      break;
      // **********************************************************************;
    case CSP_MATRIX:
      // **********************************************************************;

      break;
      // **********************************************************************;
    case ZSP_MATRIX:
      // **********************************************************************;

      break;
      // **********************************************************************;
    default:
      *ierr = blas_error_param;
      return;
    }
}
