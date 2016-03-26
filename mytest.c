#include "iBLAS.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define ROW 10
#define COL 9
#define NNZ 30

// #define ROW 6
// #define COL 4
// #define NNZ 8

// #define ROW 7
// #define COL 10
// #define NNZ 15

int main()
{
  const int N=4;
  const int nz=6;
  double val[]={1.1,2.2,2.4,3.3,4.1,4.4};
  int indx[]={0,1,1,2,3,3};
  int jndx[]={0,1,3,2,0,3};
  double x[]={1.0,1.0,1.0,1.0};
  double y[]={0.0,0.0,0.0,0.0};
  int k[4]={2,3,4,5};
  int l[3]={2,2,2};


  //for test(case) 3 lmbv_bsc
  double v3[8]={3,1,5,0,4,0,0,1};
  int bindx3[2]={0,2};
  int pb3[2]={0,1};
  int pe3[2]={1,2};
  double x3[6]={1,1,1,1,1,1};
  double y3[4];


  //for test4 rmbv_bsc
  double v4[8]={3,1,5,0,4,0,0,1};
  int bindx4[2]={0,2};
  int pb4[2]={0,1};
  int pe4[2]={1,2};
  double x4[4]={1,-1,1,-1};
  double y4[6];


  //for test5 rmbv_bsr
  double v5[8]={3,1,5,0,4,0,0,1};
  int bindx5[2]={0,1};
  int pb5[3]={0,1,1};
  int pe5[3]={1,1,2};
  double x5[4]={1,-1,1,-1};
  double y5[6];


  //for test6 lmbv_bsr
  double v6[8]={3,1,5,0,4,0,0,1};
  int bindx6[2]={0,1};
  int pb6[3]={0,1,1};
  int pe6[3]={1,1,2};
  double x6[6]={1,1,1,1,1,1};
  double y6[4];



  //for test7 lmbv_bco and rmbv_bco
  double v7[8]={3,1,5,0,4,0,0,1};
  int bindx7[2]={0,2};
  int bjndx7[2]={0,1};
  double x7[6]={1,1,1,1,1,1};
  double y7[4];
  double xx7[4]={1,-1,1,-1};
  double yy7[6];

  //for test8 lmbv_dia and rmbv_dia
  double v8[16]={0,0,4,1, 5,0,0,0, 3,0,0,0, 1,0,0,0};
  int idiag8[4]={-2,-1,0,1};
  double x8[6]={1,1,1,1,1,1};
  double y8[4];
  double xx8[4]={1,-1,1,-1};
  double yy8[6];





  //for test9 lmbv_bdi and rmbv_bdi
  double v9[16]={3,1,5,0,0,0,0,0,0,0,0,0,4,0,0,1};

  int ibdiag9[2]={0,-1};

  double x9[6]={1,1,1,1,1,1};
  double y9[4];
  double xx9[4]={1,-1,1,-1};
  double yy9[6];
  //for test10 mbv_vbr
  double v10[21]={3,1,0, 5,0,0,0,0,0, 0,0,0,0,0,4,0,0,0, 0,0,1};
  int IA1_10[4]={0,0,0,1};  
  int IA2_10[5]={0,3,9,18,21};

  int bp1_10[4]={0,1,3,6};
  int bp2_10[3]={0,3,4};

  int PB_10[3]={0,1,2};
  int PE_10[3]={1,2,4};

  double x10[6]={1,1,1,1,1,1};
  double y10[4];
  double xx10[4][1]={1,-1,1,-1};
  double yy10[6][1];


  //for test11 coo
  // double v11[NNZ]={1,5,2,6,8,3,7,9,4};
  // int IA1_11[NNZ]={0,0,1,1,2,2,2,3,3};
  // int IA2_11[NNZ]={0,1,1,2,0,2,3,1,3};
  // double v11[NNZ]={8,3,7,9,4,1,1,2};
  // int IA1_11[NNZ]={2,2,2,3,3,4,4,5};
  // int IA2_11[NNZ]={0,2,3,1,3,2,3,1};

  // double v11[NNZ]={1,1,2,8,1,1,2,3,2,3,1,4,5,9,1};
  // int IA1_11[NNZ]={1,1,3,1,6,2,2,2,1,3,3,3,4,1,6};
  // int IA2_11[NNZ]={3,0,0,8,1,1,2,3,5,2,3,7,0,9,4};

  int IA1_11[NNZ]={0,0,0,1,1,8,1,2,2,2,2,3,4,4,0,4,5,5,5,6,6,7,7,7,7,8,9,9,5,9};
  int IA2_11[NNZ]={3,4,6,0,6,3,8,2,4,5,1,6,2,4,7,5,2,3,8,5,7,0,2,4,5,5,0,2,5,8};
  double v11[NNZ]={2,9,4,1,9,4,1,2,3,4,1,2,2,3,5,5,6,4,4,8,4,2,3,4,5,3,5,2,2,8};
  
  // double xx11[4]={1,-1,1,-1};
  // double yy11[6];

  double v13[7]={3,1,1,5,1,4,1};
  int IA1_13[7]={0,0,0,1,4,4,5};
  int IA2_13[7]={0,1,3,0,0,2,3};

  double xx13[4] = {1,2,3,4};
  double yy13[6] = {1,1,1,1,1,1};

  //test level 1
  double xx12[5]={3,1,5,4,1};
  double yy12[5]={1,2,2,1,2};
  int yindex[5]={0,1,2,3,4};
  double r = 5.0;
  double re;

   // double xx10[4][1]={1,-1,1,-1};
   // double yy10[6][1];

  d_matrix* A;
  dsp_linknode* dsp_A,*dsp_B;
  DSPMAT* dsp_data;
  int i;
  double alpha=1.0;
  double begin,end;
  int result=0;
  int n = 10;
  int test = 0;

  double* x11 = (double*)malloc(sizeof(double)*ROW);
  double* x_11 = (double*)malloc(sizeof(double)*ROW);
  double* y11 = (double*)malloc(sizeof(double)*COL);
  double* y_11 = (double*)malloc(sizeof(double)*COL);
  memset(x11,0.0,sizeof(double)*ROW);
  memset(x_11,0.0,sizeof(double)*ROW);
  memset(y11,0.0,sizeof(double)*COL);
  memset(y_11,0.0,sizeof(double)*COL);
  for(i=0;i<ROW;i++){
    x11[i] = 1;
  }
  for(i=0;i<COL;i++){
    y_11[i] = 1.0;
  }

  printf("Enter the test:");

  scanf("%d",&test);

  switch(test)
    {
    case 0:
      A=BLAS_duscr_begin (N,N); 
      for(i=0;i<nz;i++)
        BLAS_duscr_insert_entry(A,val[i],indx[i],jndx[i]);
      dump_matrix (A);
      dsp_A=BLAS_uscr_end (A);
      //BLAS_dusmv (ORIGIN_MATRIX,alpha,dsp_A,x,1,y,1);
      result=BLAS_usds (dsp_A,3);
      break;
    case 1:
      A=BLAS_duscr_variable_block_begin (4,3,k,l);
      dump_matrix (A);
      begin  = omp_get_num_threads();
      for(i=0;i<nz;i++)
        BLAS_duscr_insert_entry (A,val[i],indx[i],jndx[i]);
      end = omp_get_num_threads();
      dsp_A=BLAS_uscr_end(A);
      result=BLAS_usds (dsp_A,3);
      printf("time:%.30lf\n",end - begin);
      break;
    case 2:
      result = 1;
      //test coo
      dsp_A=duscr_coo (N,N,val,nz,indx,jndx,nz,1088,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_coo(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //    printf("%lf\n",x_11[i]);

      //usconv_coo2bco(dsp_A,&result);

      //test csc
      // usconv_coo2csc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_csc(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);
      // usconv_csc2coo(dsp_A,&result);

      //test dia
      usconv_coo2dia(dsp_A,&result);
      printf("DIA format sparse matrix-vector multiplication: y = Ax\n");
      dsp_data = accessdata_dsp(dsp_A,&result);
      drmbv_dia(dsp_data,x,N,y,N,&result);
      for(i=0;i<N;i++)
        printf("%lf\n",y[i]);
      //usconv_dia2coo(dsp_A,&result);

      //test csr
      // usconv_coo2csr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_csr(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);
      // usconv_csr2coo(dsp_A,&result);   

      //test bco
      //usconv_csr2bco(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_bco(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);

      //test bsc
      // usconv_bco2bsc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_bsc(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);
      // usconv_bsc2bco(dsp_A,&result);

      //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_bsr(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);
      // usconv_bsr2bco(dsp_A,&result);

      // usconv_bco2coo(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bco(dsp_data,x11,6,y11,4,&result);
      // for(i=0;i<4;i++)
      //   printf("%lf\n",y11[i]);

      //test bdi
      // usconv_bco2bdi(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bdi(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);
      // usconv_bdi2bco(dsp_A,&result);
      
      result=BLAS_usds (dsp_A,3);
      printf("the result 2 is:%d\n",result);
      break;
    case 3:
      dsp_A=duscr_bsc (6,4,v3,8,bindx3,2,pb3,2,pe3,2,3,2,2,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_bsc (dsp_data,x3,6,y3,4,&result);
      for(i=0;i<4;i++)
        printf("%f\n",y3[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 4:
      result=1;
      dsp_A=duscr_bsc (6,4,v4,8,bindx4,2,pb4,2,pe4,2,3,2,2,1088,&result);
      usconv_bsc2bco(dsp_A,&result);
      usconv_bco2bsc(dsp_A,&result);
      drmbv_bsc(dsp_A,x4,4,y4,6,&result);

      for(i=0;i<6;i++)
        printf("%f\n",y4[i]);
      result=BLAS_usds (dsp_A,3);
      printf("the rmbv_bsc result is  %d\n",result);
      break;
    case 5:
      result=1;
      dsp_A=duscr_bsr (6,4,v5,8,bindx5,2,pb5,3,pe5,3,3,2,2,1088,&result);
      usconv_bsr2bco(dsp_A,&result);
      dsp_A->number=1;
      printf("res is %d \n",dsp_A->number);
      usconv_bco2bsr(dsp_A,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      drmbv_bsr(dsp_data,x5,4,y5,6,&result);
      for(i=0;i<6;i++)
        printf("%f\n",y5[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result 5 is  %d\n",result);
      break;
    case 6:
      dsp_A=duscr_bsr (6,4,v6,8,bindx6,2,pb6,3,pe6,3,3,2,2,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_bsr(dsp_data,x6,6,y6,4,&result);

      for(i=0;i<4;i++)
        printf("%f\n",y6[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 7:
      result=1;
      dsp_A=duscr_bco (6,4,v7,8,bindx7,2,bjndx7,2,2,3,2,2,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_bco(dsp_data,x7,6,y7,4,&result);
      for(i=0;i<4;i++)
        printf("%f\n",y7[i]);
      printf("the result is  %d\n",result);
      drmbv_bco(dsp_data,xx7,4,yy7,6,&result);
      for(i=0;i<6;i++)
        printf("%f\n",yy7[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 8:
      dsp_A=duscr_dia (6,4,v8,16,4,idiag8,4,4,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_dia(dsp_data,x8,6,y8,4,&result);
      for(i=0;i<4;i++)
        printf("%f\n",y8[i]);
      printf("the result is  %d\n",result);
      drmbv_dia(dsp_data,xx8,4,yy8,6,&result);
      for(i=0;i<6;i++)
        printf("%f\n",yy8[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 9:
      result=1;
      dsp_A=duscr_bdi (6,4,v9,16,2,ibdiag9,2,2,2,2,2,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      usconv_bdi2bco(dsp_A,&result);
      usconv_bco2bdi(dsp_A,&result);
      dlmbv_bdi(dsp_data,x9,6,y9,4,&result);
      for(i=0;i<4;i++)
        printf("%f\n",y9[i]);
      printf("the result is  %d\n",result);
      drmbv_bdi(dsp_data,xx9,4,yy9,6,&result);
      for(i=0;i<6;i++)
        printf("%f\n",yy9[i]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 10:
      dsp_A=duscr_vbr (6,4,v10,21,IA2_10,5,IA1_10,4,bp1_10,4,bp2_10,3,PB_10,3,PE_10,3,3,2,1088,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_vbr(dsp_data,x10,6,y10,4,&result);
      for(i=0;i<4;i++)
        printf("%f\n",y10[i]);
      printf("the result is~  %d\n",result);

      for(i=0;i<6;i++)
        yy10[i][0]=0;
      double* xx_10 = (double*)xx10;
      double* yy_10 = (double*)yy10;
      result=BLAS_dusmm (1,ORIGIN_MATRIX,1,1,dsp_data,xx_10,4,yy_10,6);
      for(i=0;i<6;i++)
        printf(": %f\n",yy10[i][0]);

      result=BLAS_usds (dsp_A,3);
      printf("the result is  %d\n",result);
      break;
    case 11:
      result = 1;
      //test coo
      dsp_A=duscr_coo (ROW,COL,v11,NNZ,IA1_11,IA2_11,NNZ,1088,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_coo(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);

      // usconv_coo2bco(dsp_A,&result);
      // usconv_bco2coo(dsp_A,&result);
      // usconv_coo2csr(dsp_A,&result);
      // usconv_csr2bco(dsp_A,&result);

      // //test csc
      // usconv_coo2csc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_csc(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_csc2coo(dsp_A,&result);

      // //test dia
      usconv_coo2dia(dsp_A,&result);
      printf("DIA format sparse matrix-vector multiplication\n");
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_dia(dsp_data,x11,ROW,y11,COL,&result);
      for(i=0;i<COL;i++)
        printf("%lf\n",y11[i]);
      //usconv_dia2coo(dsp_A,&result);

      //test csr
      //usconv_coo2csr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_csr(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);

      //test bco
      //usconv_csr2bco(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bco(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);

      // //test bsc
      // usconv_bco2bsc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bsc(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_bsc2bco(dsp_A,&result);

      // //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bsr(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_bsr2bco(dsp_A,&result);

      //test bdi
      // usconv_bco2bdi(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bdi(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_bdi2bco(dsp_A,&result);
      
      result=BLAS_usds (dsp_A,3);
      printf("the result 11 is:%d\n",result);
      break;
    case 12:
      result = 1;
      begin = omp_get_num_threads();
      for(i=0;i<5000000;i++){
        re = BLAS_dusdot(5,xx12,yindex,yy12,5,BLAS_BASE_ZERO);
      }
      end = omp_get_num_threads();
      printf("usdot:%lf   time:%.30lf\n",re,end - begin);
      // BLAS_dusaxpy (5,2.0,xx12,yindex,yy12,5,BLAS_BASE_ZERO);
      // printf("usaxpy:");
      // for(i=0;i<5;i++){
      //   printf("%lf\t",yy12[i]);
      // }
      
      // dsp_A = duscr_coo(6,4,v13,7,IA1_13,IA2_13,7,1088,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // BLAS_dusmm(1,ORIGIN_MATRIX,1,r,dsp_data,xx13,4,yy13,6);
      // for(i=0;i<6;i++){
      //   printf("%lf\t",yy13[i]);
      // }
      // printf("\n");
      //result = BLAS_usds(dsp_A,3);
      break;
    case 13:
      result = 1;
      dsp_A=duscr_coo (N,N,val,nz,indx,jndx,nz,1088,&result);
      usconv_coo2csr(dsp_A,&result);
      usconv_csr2coo(dsp_A,&result);
      break;

    default:
      break;
    }

  printf("hello %d\n",result);
  return 0;
}
