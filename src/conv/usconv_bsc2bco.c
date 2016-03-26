#include "properties.h"
#include "conv/conv_tools.h"
#include "link.h"
#include "conv/usconv.h"
#include <string.h>
void usconv_bsc2bco(void* spmat,int *ierr)
{
  int res,rest,i,j,kb;
  ISPMAT* isp_data;
  SSPMAT* ssp_data;
  DSPMAT* dsp_data;
  CSPMAT* csp_data;
  ZSPMAT* zsp_data;
  *ierr=-1;
  rest=((dsp_linknode*)spmat)->number;
  switch (rest){
    case ISP_MATRIX:
      // **********************************************************************

      break;
      // **********************************************************************
    case SSP_MATRIX:
      // **********************************************************************

      break;
      // **********************************************************************
    case DSP_MATRIX:
      // **********************************************************************
      dsp_data=accessdata_dsp(spmat,ierr);
      if (*ierr!=0) {
          *ierr = blas_error_param;
          return;
        }
      if(dsp_data->FIDA==BSC_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dsp_data->FIDA=BCO_FORMAT;
              get_infoa(dsp_data->INFOA,'g',&kb,ierr);
              dsp_data->IA2=(int*)aligned_malloc (sizeof(int)*(dsp_data->n_IA1));
              memset(dsp_data->IA2,0,sizeof(int)*dsp_data->n_IA1);
              dsp_data->n_IA2 = dsp_data->n_IA1;
              for(i = 0; i < kb; i++){
                  for(j = dsp_data->PB[i]; j < dsp_data->PE[i];j++){
                      dsp_data->IA2[j] = i;
                  }
              }
              if(dsp_data->PB != NULL){
                aligned_free(dsp_data->PB);
                dsp_data->PB = NULL;
                dsp_data->n_PB = 0;
              }
              if(dsp_data->PE != NULL){
                aligned_free(dsp_data->PE);
                dsp_data->PE=NULL; 
                dsp_data->n_PE = 0; 
              }

              // printf("\n--------------BSC2BCO---------------------\n");
              // printf("n_IA1:%d  n_IA2:%d  n_A:%d\n",dsp_data->n_IA1,dsp_data->n_IA2,dsp_data->n_A);
              // printf("\n");
              // for(i=0;i<dsp_data->n_IA1;i++){
              //   printf("%d\t%d",dsp_data->IA1[i],dsp_data->IA2[i]);
              // }
              // printf("\n");
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   printf("%lf\t",dsp_data->A[i]);
              // }
              // printf("\n-----------------BSC2BCO------------------\n");

            }
        }else{
          *ierr = blas_error_param;
          return;
        }
      break;
      // **********************************************************************
    case CSP_MATRIX:
      // **********************************************************************

      break;
      // **********************************************************************
    case ZSP_MATRIX:
      // **********************************************************************

      break;
      // **********************************************************************
    default:
      *ierr = blas_error_param;
      return;
    }
}//end subroutine usconv_bsc2bco
