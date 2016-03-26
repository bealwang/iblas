#define SUPPORT_MIC
#include "iBLAS.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define ROW 10
#define COL 9
#define NNZ 30

// #define ROW 7
// #define COL 10
// #define NNZ 15

double compute_time = 0.0;

int main()
{
  //for test11 coo
  // double v11[9]={1,5,2,6,8,3,7,9,4};
  // int IA1_11[9]={0,0,1,1,2,2,2,3,3};
  // int IA2_11[9]={0,1,1,2,0,2,3,1,3};
  // double v11[5]={8,3,7,9,4};
  // int IA1_11[5]={2,2,2,3,3};
  // int IA2_11[5]={0,2,3,1,3};

  // double v11[NNZ]={1,1,2,8,1,1,2,3,2,3,1,4,5,9,1};
  // int IA1_11[NNZ]={1,1,3,1,6,2,2,2,1,3,3,3,4,1,6};
  // int IA2_11[NNZ]={3,0,0,8,1,1,2,3,5,2,3,7,0,9,4};

  int IA1_11[NNZ]={0,0,0,1,1,8,1,2,2,2,2,3,4,4,0,4,5,5,5,6,6,7,7,7,7,8,9,9,5,9};
  int IA2_11[NNZ]={3,4,6,0,6,3,8,2,4,5,1,6,2,4,7,5,2,3,8,5,7,0,2,4,5,5,0,2,5,8};
  double v11[NNZ]={2,9,4,1,9,4,1,2,3,4,1,2,2,3,5,5,6,4,4,8,4,2,3,4,5,3,5,2,2,8};
  
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
    case 2:
      result = 1;
      //test coo
      dsp_A=duscr_coo (ROW,COL,v11,NNZ,IA1_11,IA2_11,NNZ,1088,&result);
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
      dsp_data = accessdata_dsp(dsp_A,&result);
      drmbv_dia(dsp_data,y_11,COL,x_11,ROW,&result);
      for(i=0;i<ROW;i++)
        printf("%lf\n",x_11[i]);
      usconv_dia2coo(dsp_A,&result);

      //test csr
      usconv_coo2csr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // drmbv_csr(dsp_data,y_11,COL,x_11,ROW,&result);
      // for(i=0;i<ROW;i++)
      //   printf("%lf\n",x_11[i]);

      //test bco
      usconv_csr2bco(dsp_A,&result);
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
      usconv_bco2bdi(dsp_A,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      drmbv_bdi(dsp_data,y_11,COL,x_11,ROW,&result);
      for(i=0;i<ROW;i++)
        printf("%lf\n",x_11[i]);
      usconv_bdi2bco(dsp_A,&result);
      
      result=BLAS_usds (dsp_A,3);
      printf("the result 2 is:%d\n",result);
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

      //test csc
      // usconv_coo2csc(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_csc(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_csc2coo(dsp_A,&result);

      //test dia
      usconv_coo2dia(dsp_A,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_dia(dsp_data,x11,ROW,y11,COL,&result);
      for(i=0;i<COL;i++)
        printf("%lf\n",y11[i]);
      usconv_dia2coo(dsp_A,&result);

      //test csr
      usconv_coo2csr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_csr(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);

      //test bco
      usconv_csr2bco(dsp_A,&result);
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

      //test bsr
      // usconv_bco2bsr(dsp_A,&result);
      // dsp_data = accessdata_dsp(dsp_A,&result);
      // dlmbv_bsr(dsp_data,x11,ROW,y11,COL,&result);
      // for(i=0;i<COL;i++)
      //   printf("%lf\n",y11[i]);
      // usconv_bsr2bco(dsp_A,&result);

      //test bdi
      usconv_bco2bdi(dsp_A,&result);
      dsp_data = accessdata_dsp(dsp_A,&result);
      dlmbv_bdi(dsp_data,x11,ROW,y11,COL,&result);
      for(i=0;i<COL;i++)
        printf("%lf\n",y11[i]);
      usconv_bdi2bco(dsp_A,&result);
      
      result=BLAS_usds (dsp_A,3);
      printf("the result 11 is:%d\n",result);
      break;

    default:
      break;
    }

  printf("hello %d\n",result);
  return 0;
}
