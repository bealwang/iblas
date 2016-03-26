#include "properties.h"
#include "conv/conv_tools.h"
#include "link.h"
#include "conv/usconv.h"
void usconv_coo2dia(void* spmat,int* ierr)
{
  int res,LDA,NDIAG,nnz,rest,i;
  int *temp;
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
      if(dsp_data->FIDA==COO_FORMAT) {
          get_infoa(dsp_data->INFOA ,'c',&res,ierr);
          if(res==COP_OF_SOURCE) {
              dpre_usconv_coo2dia ( dsp_data->M, dsp_data->K,&(dsp_data->A),dsp_data->n_IA1,&(dsp_data->IA1), dsp_data->IA2,&LDA,&NDIAG);
              if(0 == NDIAG){
                *ierr = -1;
                return;
              }
              dsp_data->FIDA=DIA_FORMAT;
              dsp_data->n_A=LDA*NDIAG;
              dsp_data->n_IA1=NDIAG;
              if(dsp_data->IA2 != NULL){
                aligned_free(dsp_data->IA2);
                dsp_data->IA2 = NULL;
                dsp_data->n_IA2 = 0;
              }
              nnz=0;
              for(i=0;i<dsp_data->n_A;i++)
                if(dsp_data->A[i]!=0)
                  nnz++;
              if(nnz<=LDA*NDIAG*0.5)
                printf("\tWarning Many zeros stored\n");
              set_infoa(dsp_data->INFOA,'n',nnz,ierr);
              set_infoa(dsp_data->INFOA,'d',LDA,ierr) ;//row-dim of val
              set_infoa(dsp_data->INFOA,'e',NDIAG,ierr); //col-dim of val

              // printf("\n-----------dia-------------------\n");
              // printf("LDA:%d  NDIAG:%d\n",LDA,NDIAG);
              // for(i=0;i<dsp_data->n_IA1;i++){
              //   printf("%d  ",dsp_data->IA1[i]);
              // }
              // printf("\n");
              // for(i=0;i<dsp_data->n_A;i++){
              //   if(i%LDA == 0){
              //     printf("\n");
              //   }
              //   printf("%lf  ",dsp_data->A[i]);
              // }
              // printf("\n--------------dia----------------\n");

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
