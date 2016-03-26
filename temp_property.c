
#include "conv/usconv.h"
#include "uscr/uscr_coo.h"
#include "blas_enum.h"
#include "link.h"
#include "comm_tools.h"
#include "sp_io.h"
#include "types.h"
#include "INSERTING.h"
#include "properties.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MV_COO
#include "host/host_lmbv_coo.h"
#include "host/host_rmbv_coo.h"
#endif

#ifdef MV_CSR
#include "host/host_lmbv_csr.h"
#include "host/host_rmbv_csr.h"
#endif

#ifdef MV_CSC
#include "host/host_lmbv_csc.h"
#include "host/host_rmbv_csc.h"
#endif

#ifdef MV_DIA
#include "host/host_lmbv_dia.h"
#include "host/host_rmbv_dia.h"
#endif

#ifdef MV_BCO
#include "host/host_lmbv_bco.h"
#include "host/host_rmbv_bco.h"
#endif

#ifdef MV_BSR
#include "host/host_lmbv_bsr.h"
#include "host/host_rmbv_bsr.h"
#endif

#ifdef MV_BSC
#include "host/host_lmbv_bsc.h"
#include "host/host_rmbv_bsc.h"
#endif

#ifdef MV_BDI
#include "host/host_lmbv_bdi.h"
#include "host/host_rmbv_bdi.h"
#endif

#include <stdio.h>
#include <assert.h>

double total = 0.0;

void Action(dsp_linknode* dsp_A,double* x,int n,double* y,int m,int *ierr){
	DSPMAT* dsp_data;
	double begin,end;
	int i;
	#ifdef MV_COO
		dsp_data = accessdata_dsp(dsp_A,ierr);
		begin = omp_get_wtime();
		dlmbv_coo(dsp_data, x, n, y, m, ierr);
		end = omp_get_wtime();
	#endif
    #ifdef MV_CSR
        usconv_coo2csr(dsp_A,ierr);
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_csr(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_CSC
        usconv_coo2csc(dsp_A,ierr);
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_csc(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_DIA
        usconv_coo2dia(dsp_A,ierr);
        if(*ierr == -1){
        printf("too many diags!\n");
        exit(1);
        }
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_dia(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_BCO
        usconv_coo2csr(dsp_A,ierr);
        usconv_csr2bco(dsp_A,ierr);
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_bco(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_BSR
        usconv_coo2csr(dsp_A,ierr);
        usconv_csr2bco(dsp_A,ierr);
        usconv_bco2bsr(dsp_A,ierr);
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_bsr(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_BSC
        usconv_coo2csr(dsp_A,ierr);
        usconv_csr2bco(dsp_A,ierr);
        usconv_bco2bsc(dsp_A,ierr);
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_bsc(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
    #ifdef MV_BDI
        usconv_coo2csr(dsp_A,ierr);
        usconv_csr2bco(dsp_A,ierr);
        usconv_bco2bdi(dsp_A,ierr);
        if(*ierr == -1){
            printf("too many diags!\n");
            exit(1);
        }
        dsp_data = accessdata_dsp(dsp_A,ierr);
        begin = omp_get_wtime();
        dlmbv_bdi(dsp_data, x, n, y, m,ierr);
        end = omp_get_wtime();
    #endif
	total += end - begin;
}

int main(int argc, char *argv[]){
    int M, N, nz;   
    int i;
    char* file_name;
    int result = 1;
    dsp_linknode* dsp_A;
  	DSPMAT* dsp_data;

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
	    nz = dsp_data->n_A;
	  }

	double* x11 = (double*)malloc(sizeof(double)*M);
	double* y11 = (double*)malloc(sizeof(double)*N);
	memset(x11,0.0,sizeof(double)*M);
	memset(y11,0.0,sizeof(double)*N);

	for(i=0;i<M;i++){
	  x11[i] = 0.1;
	}
    
	for(i=0;i<500;i++){
		Action(dsp_A,x11,M,y11,N,&result);
	}

	FILE* fp = NULL;
	fp = fopen("./timer","w");
	if(fp == NULL){
		printf("error!\n");
	}
	fprintf(fp,"time:%lf\n",total);
	fclose(fp);
	fp = NULL;

	// for(i=0;i<COL;i++)
 //        printf("%lf\n",y11[i]);
    printf("hello %d\n",result);
    return 0;
}