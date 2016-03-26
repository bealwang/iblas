#include "properties.h"
#include "conv/conv_tools.h"
#include "link.h"
#include "conv/usconv.h"
#include <string.h>
//if call this funciton ,the res in dsp_data->INFOA must be COP_OF_SOURSE
void usconv_bco2bsr(void* spmat,int* ierr)
{
  int res,mb,kb,lb,bnnz,rest,i;
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
      /*INT TYPE convert*/
      break;
      // **********************************************************************;
    case SSP_MATRIX:
      // **********************************************************************;
      /*DOUBLE FLOAT TYPE convert*/
      break;
      // **********************************************************************;
    case DSP_MATRIX:
      // **********************************************************************;

      dsp_data=accessdata_dsp(spmat,ierr);
      if (*ierr!=0) {
          *ierr = blas_error_param;
          return;
        }
      if(dsp_data->FIDA==BCO_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dsp_data->FIDA=BSR_FORMAT;
              get_infoa(dsp_data->INFOA,'n',&bnnz,ierr);
              get_infoa(dsp_data->INFOA,'f',&mb,ierr);
              get_infoa(dsp_data->INFOA,'g',&kb,ierr);
              get_infoa(dsp_data->INFOA,'e',&lb,ierr);
              dsp_data->PB=(int* )aligned_malloc(sizeof(int)*mb);
              dsp_data->PE=(int* )aligned_malloc(sizeof(int)*mb);
              memset(dsp_data->PB,0,sizeof(int)*mb);
              memset(dsp_data->PE,0,sizeof(int)*mb);
              dsp_data->n_PB = mb;
              dsp_data->n_PE = mb;
              dpre_usconv_bco2bsr(&(dsp_data->IA1),dsp_data->IA2,&(dsp_data->A),mb,kb,bnnz,lb,dsp_data->PB,dsp_data->PE);
              if(dsp_data->IA2 != NULL){
                aligned_free(dsp_data->IA2);
                dsp_data->IA2 = NULL;
                dsp_data->n_IA2 = 0;
              }
              
              // printf("\n---------------------BCO2BSR-----------------------\n");
              // for(i=0;i<dsp_data->n_PE;i++){
              //   printf("PB:%d  PE:%d\n",dsp_data->PB[i],dsp_data->PE[i]);
              // }
              // for(i=0;i<dsp_data->n_IA1;i++){
              //   printf("%d\t",dsp_data->IA1[i]);
              // }
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   printf("%lf  ",dsp_data->A[i]);
              // }
              // printf("\n---------------------BCO2BSR-----------------------\n");
            }
        }else{
          *ierr = blas_error_param;
          return;
        }
      break;
      // **********************************************************************;
    case CSP_MATRIX:
      // **********************************************************************;
      /*COMPLEX SINGLE FLOAT TYPE convert*/
      break;
      // **********************************************************************;
    case ZSP_MATRIX:
      // **********************************************************************;
      /*COMPLEX DOUBLE FLOAT TYPE convert*/
      break;
      // **********************************************************************;
    default:
      *ierr = blas_error_param;
      return;
    }
}
