#include "properties.h"
#include "conv/conv_tools.h"
#include "conv/usconv.h"
#include "link.h"
void usconv_bco2bdi(void* spmat,int* ierr)
{
  int  res,BLDA,BNDIAG,mb,kb,lb,rest,i;
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
      if(dsp_data->FIDA==BCO_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              get_infoa(dsp_data->INFOA ,'e',&lb,ierr);
              get_infoa(dsp_data->INFOA ,'f',&mb,ierr);
              get_infoa(dsp_data->INFOA ,'g',&kb,ierr);
              BLDA=mb<kb?mb:kb;
              dpre_usconv_bco2bdi(mb,kb,lb,&(dsp_data->A),&(dsp_data->IA1),dsp_data->n_IA1,dsp_data->IA2,&BLDA,&BNDIAG);
              if(0 == BNDIAG){
                *ierr = -1;
                return;
              }
              dsp_data->FIDA=BDI_FORMAT;
              dsp_data->n_IA1 = BNDIAG;
              if(dsp_data->IA2 != NULL){
                aligned_free(dsp_data->IA2);
                dsp_data->IA2 = NULL;
                dsp_data->n_IA2 = 0;
              }
              dsp_data->n_A = BLDA*BNDIAG*lb*lb;
              set_infoa(dsp_data->INFOA,'d',lb,ierr); //row-dim of a block;
              set_infoa(dsp_data->INFOA,'e',lb,ierr); //col-dim of a block;
              set_infoa(dsp_data->INFOA,'f',BLDA,ierr); //blocks per diagonal;
              set_infoa(dsp_data->INFOA,'g',BNDIAG,ierr); //no of diagonals;

              // printf("\n--------------------bdi------------------------\n");
              // printf("lb:%d\tBLDA:%d\tBNDIAG:%d\n",lb,BLDA,BNDIAG);
              // for(i=0;i<dsp_data->n_IA1;i++){
              //   printf("%d\t",dsp_data->IA1[i]);
              // }
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   printf("%lf\t",dsp_data->A[i]);
              // }
              // printf("\n--------------------bdi------------------------\n");
              
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
