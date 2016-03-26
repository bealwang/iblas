#include "iBLAS.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NUMBER 500

int main(int argc, char *argv[])
{
  int i;
  double compute_time,compute_begin,compute_end,begin,end;
  char* file_name;
  double v11[5]={3,1,5,4,1};
  int IA1_11[5]={0,0,1,4,5};
  int IA2_11[5]={0,1,0,2,3};
  double xx11[4]={1,-1,1,-1};
  double yy11[6];
  //test level 1
  double xx12[5]={3,1,5,4,1};
  double yy12[5]={1,2,2,1,2};
  int yindex[5]={0,1,2,3,4};
  double r = 5.0;
  double re;

  //  double xx10[4][1]={1,-1,1,-1};
  //  double yy10[6][1];

  d_matrix* A;
  dsp_linknode* dsp_A;
  DSPMAT* dsp_data;
  double alpha=1.0;
  int result=1;
  int test,nnz,M,N;

  if (argc < 2){
    fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
    exit(1);
  }else{ 
    file_name = argv[1]; 
    dsp_A = read_coo_matrix(file_name);
  }
  if(NULL != dsp_A){
    dsp_data = accessdata_dsp(dsp_A,&result);
  }
  if(0 == result){
    M = dsp_data->M;
    N = dsp_data->K;
    nnz = dsp_data->n_A;
  }

  double* x11 = (double*)malloc(sizeof(double)*M);
  double* x_11 = (double*)malloc(sizeof(double)*M); 
  double* y11 = (double*)malloc(sizeof(double)*N);
  double* y_11 = (double*)malloc(sizeof(double)*N);
  memset(x11,0.0,sizeof(double)*M);
  memset(y11,0.0,sizeof(double)*N);
  memset(x_11,0.0,sizeof(double)*M);
  memset(y_11,0.0,sizeof(double)*N);

  for(i=0;i<M;i++){
    x11[i] = 0.1;
  }
  for(i=0;i<N;i++){
    y_11[i] = 0.1;
  }
  printf("%d  %d\n",M,N);
  printf("Enter the test:");
  scanf("%d",&test);

  switch(test)
    {
    case 10:
      result=1;
      //test coo
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_coo(dsp_data,y_11,N,x_11,M,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = (double)(end - begin);
      // conv_time = (double)(conv_end - conv_begin);
      // printf("coo_format\ncompute time:%.30lf    conv_time:%.30lf\n",compute_time,conv_time);

      //test csc
      // usconv_coo2csc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // compute_begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_csc(dsp_data,y_11,N,x_11,M,&result);
      // }
      // compute_end = omp_get_wtime();
      // compute_time = compute_end - compute_begin;
      // usconv_csc2coo(dsp_A,&result);
      // printf("csc_format\ncompute time:%.30lf\n",compute_time);

      //test csr
      // usconv_coo2csr (dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_csr(dsp_data,y_11,N,x_11,M,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = end - begin;
      // printf("csr_format\ncompute time:%.30lf\n",compute_time);

      //test bco
      // usconv_csr2bco(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_bco(dsp_data,y_11,N,x_11,M,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = end - begin;
      // printf("bco_format\ncompute time:%.30lf\n",compute_time);

      //test bsc
      // usconv_bco2bsc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_bsc(dsp_data,y_11,N,x_11,M,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = (double)(end - begin);
      // usconv_bsc2bco(dsp_A,&result);
      // printf("bsc_format\ncompute time:%.30lf\n",compute_time);

      //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   drmbv_bsr(dsp_data,y_11,N,x_11,M,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = (double)(end - begin);
      // usconv_bsr2bco(dsp_A,&result);
      // printf("bsr_format\ncompute time:%.30lf\n",compute_time);

      break;

    case 11:
      result=1;
      //test coo
      // compute_begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_coo(dsp_data,x11,M,y11,N,&result);
      // }
      // compute_end = omp_get_wtime();
      // compute_time = compute_end - compute_begin;
      // printf("coo_format\ncompute time:%.16g \n",compute_time);
    
      //test csc
      // usconv_coo2csc (dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // compute_begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_csc(dsp_data,x11,M,y11,N,&result);
      // }
      // compute_end = omp_get_wtime();
      // compute_time = compute_end - compute_begin;
      // usconv_csc2coo(dsp_A,&result);
      // printf("csc_format\ncompute time:%.16g\n",compute_time);
    
      //test dia
      usconv_coo2dia (dsp_A,&result);
      if(result != -1){
        dsp_data = accessdata_dsp(dsp_A,&result);
        compute_begin = omp_get_wtime();
        for(i=0;i<NUMBER;i++){
          dlmbv_dia(dsp_data,x11,M,y11,N,&result);
        }
        compute_end = omp_get_wtime();
        compute_time = compute_end - compute_begin;
        usconv_dia2coo(dsp_A,&result);
        printf("dia_format\ncompute time:%lf\n",compute_time);
      }
    
      //test csr
      // usconv_coo2csr (dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // compute_begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_csr(dsp_data,x11,M,y11,N,&result);
      // }
      // compute_end = omp_get_wtime();
      // compute_time = compute_end - compute_begin;
      // printf("csr_format\ncompute time:%.16g\n",compute_time);
    
      // //test bco
      // usconv_csr2bco(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      //  for(i=0;i<NUMBER;i++){
      //   dlmbv_bco(dsp_data,x11,M,y11,N,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = end - begin;
      // printf("bco_format\ncompute time:%.30lf\n",compute_time);
    
      // //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_bsr(dsp_data,x11,M,y11,N,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = end - begin;
      // usconv_bsr2bco(dsp_A,&result);
      // printf("bsr_format\ncompute time:%.30lf\n",compute_time);

      // //test bsc
      // usconv_bco2bsc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // begin = omp_get_wtime();
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_bsc(dsp_data,x11,M,y11,N,&result);
      // }
      // end = omp_get_wtime();
      // compute_time = end - begin;
      // usconv_bsc2bco(dsp_A,&result);
      // printf("bsc_format\ncompute time:%lf\n",compute_time);

      // //test bdi
      // usconv_bco2bdi(dsp_A,&result);
      // if(result != -1){
      //   dsp_data = accessdata_dsp(dsp_A,&result);
      //   begin = omp_get_wtime();
      //   for(i=0;i<NUMBER;i++){
      //     dlmbv_bdi(dsp_data,x11,M,y11,N,&result);
      //   }
      //   end = omp_get_wtime();
      //   compute_time = end - begin;
      //   usconv_bdi2bco(dsp_A,&result);
      //   printf("bdi_format\ncompute time:%.30lf\n",compute_time);
      // }

      free(x11);
      free(x_11);
      free(y_11);
      free(y11);
      break;

    case 12:
      re = BLAS_dusdot(5,xx12,yindex,yy12,5,BLAS_BASE_ZERO);
      printf("usdot:%lf\n",re);
      break;

    default:
      break;
    }
  printf("hello %d\n",result);
  return 0;
}