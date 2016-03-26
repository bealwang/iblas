#include "properties.h"
#include "conv/conv_tools.h"
#include "conv/usconv.h"
#include "link.h"
#include "malloc.h"
#include <string.h>
void usconv_coo2csr(void* spmat,int* ierr)
{
  int res,rest,i,j,nnz;
  int*  temp;
  ISPMAT*  isp_data;
  SSPMAT*  ssp_data;
  DSPMAT*  dsp_data;
  CSPMAT*  csp_data;
  ZSPMAT*  zsp_data;
  (*ierr)=-1;

  rest=((dsp_linknode*)spmat)->number;
  switch (rest){
    case ISP_MATRIX:

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
      if(dsp_data->FIDA==COO_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dsp_data->FIDA=CSR_FORMAT;
              get_infoa(dsp_data->INFOA ,'n',&nnz,ierr);
              dsp_data->PB=(int*)aligned_malloc(sizeof(int)*(dsp_data->M));     //K-->M  modified at 2015/1/9 by wang
              dsp_data->PE=(int*)aligned_malloc(sizeof(int)*dsp_data->M);
              memset(dsp_data->PB,0,sizeof(int)*(dsp_data->M));
              memset(dsp_data->PE,0,sizeof(int)*(dsp_data->M));
              dsp_data->n_PB=dsp_data->M;
              dsp_data->n_PE=dsp_data->M;
              //dpre_usconv_coo2csr(dsp_data->A,dsp_data->IA1,dsp_data->n_IA1,dsp_data->IA2,dsp_data->M,dsp_data->PB,dsp_data->PE);
              dpre_usconv_coo2csr(&(dsp_data->IA1),dsp_data->IA2,&(dsp_data->A),dsp_data->M,dsp_data->K,nnz,dsp_data->PB,dsp_data->PE);
              if(dsp_data->IA2 != NULL){
                aligned_free(dsp_data->IA2);
                dsp_data->IA2 = NULL;
                dsp_data->n_IA2=0;
              }

              printf("\n-------------------COO2CSR:----------------\n");
              for(i=0;i<dsp_data->n_PB;i++){
                printf("PB:%d  PE:%d\n",dsp_data->PB[i],dsp_data->PE[i]);
              }
              printf("cols:");
              for(i=0;i<dsp_data->n_IA1;i++){
                printf("%d  ",dsp_data->IA1[i]);
              }
              printf("\n");
              for(i=0;i<dsp_data->n_PB;i++){
                for(j=dsp_data->PB[i];j<dsp_data->PE[i];j++){
                  printf("%lf  ",dsp_data->A[j]);
                }
                printf("\n");
              }
              printf("\n-------------------COO2CSR:----------------\n");
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
