#include "properties.h"
#include "conv/conv_tools.h"
#include "link.h"
#include "conv/usconv.h"
#include <string.h>
void usconv_csr2coo(void* spmat,int* ierr)
{
  int res,rest,i,j;
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
      if(dsp_data->FIDA==CSR_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dsp_data->FIDA=COO_FORMAT;
              dsp_data->IA2=(int*)aligned_malloc(sizeof(int)*dsp_data->n_A);
              memset(dsp_data->IA2,0,sizeof(int)*dsp_data->n_A);
              dsp_data->n_IA2=dsp_data->n_A;
              for(i=0;i<dsp_data->n_IA1;i++)
                (dsp_data->IA2)[i]=(dsp_data->IA1)[i];
              for(i = 0; i < dsp_data->M; i++){
                  for(j = dsp_data->PB[i]; j < dsp_data->PE[i];j++){
                      dsp_data->IA1[j] = i;
                  }
              }
              dsp_data->n_BP1=0;
              dsp_data->n_BP2=0;
              if(dsp_data->PB != NULL){
                aligned_free(dsp_data->PB);
                dsp_data->PB = NULL;
                dsp_data->n_PB = 0;
              }
              if(dsp_data->PE != NULL){
                aligned_free(dsp_data->PE);
                dsp_data->PE = NULL;
                dsp_data->n_PE = 0;
              }
              

              printf("\n-----------------CSR2COO:----------------\n");
              for(i=0;i<dsp_data->n_IA1;i++){
                printf("row:%d  col:%d  value:%lf\n",dsp_data->IA1[i],dsp_data->IA2[i],dsp_data->A[i]);
              }
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   printf("%lf  ",dsp_data->A[i]);
              // }
              printf("\n----------------CSR2COO:----------------\n");

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
