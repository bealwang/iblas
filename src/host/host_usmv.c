#include "host/host_USMV.h"
#include "malloc.h"
int BLAS_dusmv (blas_trans_type transa, double alpha, dsp_linknode *A, double *x, int incx, double *y, int incy)
{
  int ierr;
  double* z;
  DSPMAT* dspmtx;
  int transa_work;
  double alpha_work;
  int i;
  ierr=-1;
  transa_work = transa;

  //defalut value should: transa_work = ORIGIN_MATRIX
  alpha_work=alpha;
  if (alpha_work==0)
    printf("no matrix multiplication necessary-------said by functin:  usmv \n");
  else
    {
      dspmtx=accessdata_dsp(A,&ierr);
      if (ierr!=0)
        {
          ierr = blas_error_param;
          return ierr;;
        }
      z=(double*)aligned_malloc(sizeof(double)*incy);

      if (ierr!=0)
        {
          ierr = blas_error_memalloc;
          return ierr;
        }

      if(transa_work==ORIGIN_MATRIX){
          switch (dspmtx->FIDA){
            case COO_FORMAT:
              drmbv_coo(dspmtx,x,incx,z,incy,&ierr);
              break;
            case CSC_FORMAT:
              drmbv_csc(dspmtx,x,incx,z,incy,&ierr);
              break;
            case CSR_FORMAT:
              drmbv_csr(dspmtx,x,incx,z,incy,&ierr);
              break;
            case DIA_FORMAT:
              drmbv_dia(dspmtx,x,incx,z,incy,&ierr);
              break;
            case BCO_FORMAT:
              drmbv_bco(dspmtx,x,incx,z,incy,&ierr);
              break;
            case BSC_FORMAT:
              drmbv_bsc(dspmtx,x,incx,z,incy,&ierr);
              break;
            case BSR_FORMAT:
              drmbv_bsr(dspmtx,x,incx,z,incy,&ierr);
              break;
            case BDI_FORMAT:
              drmbv_bdi(dspmtx,x,incx,z,incy,&ierr);
              break;
            case VBR_FORMAT:
              drmbv_vbr(dspmtx,x,incx,z,incy,&ierr);
              break;
            default:
              ierr = blas_error_param;
            }
        }
          else if(transa_work==TRANSP_MATRIX)
            {
              switch (dspmtx->FIDA){
                case COO_FORMAT:
                  dlmbv_coo(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case CSC_FORMAT:
                  dlmbv_csc(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case CSR_FORMAT:
                  dlmbv_csr(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case DIA_FORMAT:
                  dlmbv_dia(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case BCO_FORMAT:
                  dlmbv_bco(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case BSC_FORMAT:
                  dlmbv_bsc(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case BSR_FORMAT:
                  dlmbv_bsr(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case BDI_FORMAT:
                  dlmbv_bdi(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                case VBR_FORMAT:
                  dlmbv_vbr(dspmtx,(x),incx,z,incy,&ierr);
                  break;
                default:
                  ierr = blas_error_param;
                  break;
                }
        }
              else
              ierr = blas_error_param;

              if (ierr!=0)
                {
                  aligned_free(z);
                  return ierr;
                }
              if(transa_work==ORIGIN_MATRIX)
                for(i=0;i<incy;i++)
                  y[i]= alpha_work * z[i] + y[i];
              else
                for(i=0;i<incy;i++)
                  y[i]= alpha_work * ( (z[i])) + y[i];
              aligned_free(z);
            }
          ierr = 0;
          return ierr;
        }
