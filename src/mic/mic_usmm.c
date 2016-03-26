#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mic/mic_USMM.h"
int BLAS_dusmm (blas_order_type order, blas_trans_type transa, int nrhs, double alpha, DSPMAT* A, double* B, int ldb, double*  C, int ldc)
{

  DSPMAT* dspmtx;
  int transa_work;
  double *bb,*cc;
  int i,j;
  double alpha_work;
  int ierr = -1;

  alpha_work=alpha;
  transa_work=transa;
  if (alpha_work==0)
    printf("no matrix multiplication necessary-------said by functin:  usmv \n");
  else
    {
      dspmtx=A;    //dspmtx=accessdata_dsp(A,&ierr) 修改
      ierr = 0;
      if (ierr!=0)
        {
          ierr = blas_error_param;
          return ierr;;
        }
      bb=(double*)aligned_malloc(sizeof(double)*ldb);
      cc=(double*)aligned_malloc(sizeof(double)*ldc);
      if (!bb||!cc)
        {
          ierr = blas_error_memalloc;
          return ierr;;
        }
      for(i=0;i<nrhs;i++)
        {
          for(j=0;j<ldb;j++)
            bb[j]=B[j*nrhs+i];
          if(transa_work==ORIGIN_MATRIX)
            {
              switch (dspmtx->FIDA){
                case COO_FORMAT:
                  drmbv_coo(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case CSC_FORMAT:
                  drmbv_csc(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case CSR_FORMAT:
                  drmbv_csr(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case DIA_FORMAT:
                  drmbv_dia(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case BCO_FORMAT:
                  drmbv_bco(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case BSC_FORMAT:
                  drmbv_bsc(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case BSR_FORMAT:
                  drmbv_bsr(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case BDI_FORMAT:
                  drmbv_bdi(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                case VBR_FORMAT:
                  drmbv_vbr(dspmtx,bb,ldb,cc,ldc,&ierr);
                  break;
                default:
                  ierr = blas_error_param;
                }
            }
          else if(transa_work==TRANSP_MATRIX)
            {
              switch (dspmtx->FIDA){
                case COO_FORMAT:
                  dlmbv_coo(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case CSC_FORMAT:
                  dlmbv_csc(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case CSR_FORMAT:
                  dlmbv_csr(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case DIA_FORMAT:
                  dlmbv_dia(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case BCO_FORMAT:
                  dlmbv_bco(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case BSC_FORMAT:
                  dlmbv_bsc(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case BSR_FORMAT:
                  dlmbv_bsr(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case BDI_FORMAT:
                  dlmbv_bdi(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                case VBR_FORMAT:
                  dlmbv_vbr(dspmtx,(bb),ldb,cc,ldc,&ierr);
                  break;
                default:
                  ierr = blas_error_param;
                  break;
                }
            }
          else
            ierr = blas_error_param;


          if(ierr!=0){
            aligned_free(bb);
            aligned_free(cc);
            return ierr;
          }

          if(transa_work==ORIGIN_MATRIX)
            for(j=0;j<ldc;j++)
              C[j*nrhs+i]+=alpha_work*cc[j];
          else
            for(j=0;j<ldc;j++)
              C[j*nrhs+i]+=alpha_work*(cc[j]);
        }
      aligned_free(bb);
      aligned_free(cc);
    }


  ierr = 0;
  return ierr;

}
