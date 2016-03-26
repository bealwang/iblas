#include "iBLAS.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "sp_io.h"

#define NUMBER 500

double compute_time = 0.0; 

int main(int argc, char *argv[])
{
  int i;
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
  printf("Enter the test:");
  scanf("%d",&test);

  switch(test)
    {
    case 10:
      // usconv_coo2dia (dsp_A,&result);
      // if(result != -1){
      //   dsp_data = accessdata_dsp(dsp_A,&result);
      //   for(i=0;i<NUMBER;i++){
      //     drmbv_dia(dsp_data,y_11,N,x_11,M,&result);
      //   }
      //   usconv_dia2coo(dsp_A,&result);
      //   printf("dia_format\ncompute time:%lf\n",compute_time);
      // }
      usconv_coo2csr (dsp_A,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      for(i=0;i<NUMBER;i++){
        drmbv_csr(dsp_data,y_11,N,x_11,M,&result);
      }
      printf("csr_format\ncompute time:%.30lf\n",compute_time);
    
      
      // usconv_csr2bco(dsp_A,&result);

      // //test bdi
      // usconv_bco2bdi(dsp_A,&result);
      // if(result != -1){
      //   dsp_data = accessdata_dsp(dsp_A,&result);
      //   for(i=0;i<NUMBER;i++){
      //     drmbv_bdi(dsp_data,y_11,N,x_11,M,&result);
      //   }
      //   usconv_bdi2bco(dsp_A,&result);
      //   printf("bdi_format\ncompute time:%.30lf\n",compute_time);
      // }

      free(x11);
      free(x_11);
      free(y_11);
      free(y11);

      break;
    case 11:
      result=1;
      //test coo
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_coo(dsp_data,x11,M,y11,N,&result);
      // }
      // printf("coo_format\ncompute time:%.16g \n",compute_time);
    
      //test csc
      // usconv_coo2csc (dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_csc(dsp_data,x11,M,y11,N,&result);
      // }
      // usconv_csc2coo(dsp_A,&result);
      // printf("csc_format\ncompute time:%.30lf\n",compute_time);
    
      //test dia
      usconv_coo2dia (dsp_A,&result);
      if(result != -1){
        dsp_data = accessdata_dsp(dsp_A,&result);
        for(i=0;i<NUMBER;i++){
          dlmbv_dia(dsp_data,x11,M,y11,N,&result);
        }
        usconv_dia2coo(dsp_A,&result);
        printf("dia_format\ncompute time:%lf\n",compute_time);
      }
    
      //test csr
      // usconv_coo2csr (dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_csr(dsp_data,x11,M,y11,N,&result);
      // }
      // printf("csr_format\ncompute time:%.30lf\n",compute_time);
    
      //test bco
      //usconv_csr2bco(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      //  for(i=0;i<NUMBER;i++){
      //   dlmbv_bco(dsp_data,x11,M,y11,N,&result);
      // }
      // printf("bco_format\ncompute time:%.30lf\n",compute_time);
    
      //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_bsr(dsp_data,x11,M,y11,N,&result);
      // }
      // usconv_bsr2bco(dsp_A,&result);
      // printf("bsr_format\ncompute time:%.30lf\n",compute_time);

      //test bsc
      // usconv_bco2bsc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // for(i=0;i<NUMBER;i++){
      //   dlmbv_bsc(dsp_data,x11,M,y11,N,&result);
      // }
      // usconv_bsc2bco(dsp_A,&result);
      // printf("bsc_format\ncompute time:%.30lf\n",compute_time);

      // //test bdi
      // usconv_bco2bdi(dsp_A,&result);
      // if(result != -1){
      //   dsp_data = accessdata_dsp(dsp_A,&result);
      //   for(i=0;i<NUMBER;i++){
      //     dlmbv_bdi(dsp_data,x11,M,y11,N,&result);
      //   }
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