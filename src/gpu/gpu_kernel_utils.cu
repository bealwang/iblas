#include <stdio.h>
#include <time.h>
#include <cuda_runtime.h>
#include <assert.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include "properties.h"

extern double compute_time;

void* malloc_gpu(const size_t size){
	void* ptr;
	if(size > 0){
		cudaMalloc((void**)&ptr,size);
	}
	assert(ptr != NULL);
	return ptr;
}

void free_gpu(void* ptr){
	assert(ptr != NULL);
	cudaFree(ptr);
}

__global__ void mbv_coo(int* row,int* col,double* val,int nnz,double* x,double* y){
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	while(tid < nnz){
		y[col[tid]] += val[tid]*x[row[tid]];
		__syncthreads();
		tid += blockDim.x * gridDim.x;
	}
}

__global__ void mbv_csc(int* row,double* val,int* PB,int* PE,int ncol,double* x,double* y,char tag){
	int index = blockDim.x*blockIdx.x+threadIdx.x;
	if('l' == tag){
		if(index < ncol){
			int col_begin = PB[index];
			int col_end = PE[index];
			int i;
			double sum = 0.0;
			for (i=col_begin;i<col_end;i++)
			{
			  sum+=val[i]*x[row[i]];
			}
			y[index] = sum;
		}
	}else{
		if(index<ncol){
			int col_begin = PB[index];
			int col_end = PE[index];
			int i;
			for(i=col_begin;i<col_end;i++){
			  y[row[i]]+=val[i]*x[index];
			}
        }
	}
}

__global__ void mbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y,char tag){
	int index = blockDim.x*blockIdx.x+threadIdx.x;
	if(tag == 'l'){
		if(index < ncol){
			double sum = 0.0;
			int j;
			for (j=0; j<ndiag; ++j){
				int offset = diag[j];
				int istart = offset >= 0?0:(-offset);
				int jstart = offset >= 0?offset:0;
				int N = 0; int v_offset = 0; int I_offset = 0;

				if((nrow-istart) > (ncol-jstart)){
					N = ncol - jstart;
					I_offset = istart - jstart;
					v_offset = j*lda - jstart;
				}else{
					N = nrow - istart;
					I_offset = istart - jstart;
					v_offset = j*lda + lda - jstart - N;
				}

				int jend = jstart + N;

				if ((index >= jstart) && (index < jend)){
					sum += val[index+v_offset] * x[index+I_offset];
				}
				if(index >= jend){
					continue;
				}
			}
			y[index] = sum;
		}
	}else{
		if(index < nrow){
          double sum = 0.0;
          int j;
          for (j=0; j<ndiag; ++j){
            int offset = diag[j];
            int istart = offset >= 0?0:(-offset);
            int jstart = offset >= 0?offset:0;
            int N = 0; int v_offset = 0; int J_offset = 0;

            if((nrow-istart) >= (ncol-jstart)){
              N = ncol - jstart;
              J_offset = jstart - istart;
              v_offset = j*lda - istart;
            }else{
              N = nrow - istart;
              J_offset = jstart - istart;
              v_offset = j*lda + lda - nrow;
            }

            int iend = istart + N;
            
            if ((index >= istart) && (index < iend)){
              sum += val[index+v_offset] * x[index+J_offset];
            }
            if(index >= iend){
              continue;
            }
          }
          y[index] = sum;
        }
	}
}

__global__ void mbv_bsc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nb,int lda,double* x,double*y){
	int cols = blockDim.x*blockIdx.x+threadIdx.x;
	int nn = lda;
    if(cols < nb){
    	int col_begin = PB[cols];
		int col_end = PE[cols];
		int i;
		for(i=col_begin;i<col_end;i++)
		{
		  int p;
		  for (p=0;p<nn;p++){
		    int q;
		    for (q=0;q<nn;q++){
		        if((cols*nn+q < ncol) && (row[i]*nn+p < nrow)){
		          y[cols*nn+q]+=val[i*nn*nn+p*nn+q]*x[row[i]*nn+p];
		        }
		    }
		  }
		}
    }        
}

__global__ void mbv_bdi(int*bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y,char tag){
	int index = blockDim.x*blockIdx.x+threadIdx.x;
	if('l' == tag){
		if(index<nb) {
			int j;
			for (j=0; j<nbdiag; ++j){
			  int offset = bdiag[j];
			  int istart = offset >= 0?0:(-offset);
			  int jstart = offset >= 0?offset:0;
			  int N = 0; int v_offset = 0; int I_offset = 0;

			  if((mb-istart) > (nb-jstart)){
			    N = nb - jstart;
			    I_offset = istart - jstart;
			    v_offset = j*blda - jstart;
			  }else{
			    N = mb - istart;
			    I_offset = istart - jstart;
			    v_offset = j*blda + blda - jstart - N;
			  }

			  int jend = jstart + N;
			  
			  if ((index >= jstart) && (index < jend)){
			    int p;
			    for (p=0;p<nn;p++){
			      int q;
			      for (q=0;q<nn;q++){
			        if((index*nn+q < ncol) && ((index+I_offset)*mm+p < nrow)){
			          y[index*nn+q] += val[(index+v_offset)*mm*nn+nn*p+q] * x[(index+I_offset)*mm+p];
			        }
			      }
			    }
			  }
			  if(index >= jend){
			    continue;
			  }
			}
		}
	}else{
		if(index<mb) {
			int j;
			for (j=0; j<nbdiag; ++j){
			  int offset = bdiag[j];
			  int istart = offset >= 0?0:(-offset);
			  int jstart = offset >= 0?offset:0;
			  int N = 0; int v_offset = 0; int J_offset = 0;

			  if((mb-istart) >= (nb-jstart)){
			    N = nb - jstart;
			    J_offset = jstart - istart;
			    v_offset = j*blda - istart;
			  }else{
			    N = mb - istart;
			    J_offset = jstart - istart;
			    v_offset = j*blda + blda - mb;
			  }

			  int iend = istart + N;
			  
			  if ((index >= istart) && (index < iend)){
			    int p;
			    for (p=0;p<mm;p++){
			      int q;
			      for (q=0;q<nn;q++){
			        if((index*mm+p <nrow) && ((index+J_offset)*nn+q < ncol)){
			          y[index*mm+p] += val[(index+v_offset)*mm*nn+nn*p+q] * x[(index+J_offset)*nn+q];
			        }
			      }
			    }
			  }
			  if(index >= iend){
			    continue;
			  }
			}
		}
	}
}

__global__ void mbv_csr(int* row,double* val,int* PB,int* PE,int nrow,double* x,double* y){
	int index = blockDim.x*blockIdx.x+threadIdx.x;
	if(index < nrow){
		int row_begin = PB[index];
		int row_end = PE[index];
		int i;
		double sum = 0.0;
		for (i=row_begin;i<row_end;i++)
		{
		  y[index]+=val[i]*x[row[i]];
		}
	}
}

extern "C" void lmbv_coo(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y){
	int *dev_row,*dev_col;
	double *dev_val,*dev_x,*dev_y;

	dev_row = (int*)malloc_gpu(nnz*sizeof(int));
	dev_col = (int*)malloc_gpu(nnz*sizeof(int));
	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_row,row,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_col,col,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);

	mbv_coo<<<(nnz+4-1)/4,4>>>(dev_row,dev_col,dev_val,nnz,dev_x,dev_y);

	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row);
	free_gpu(dev_col);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void lmbv_csr(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y){
	cusparseHandle_t cusparseHandle = 0;
	cusparseMatDescr_t mat_descr = 0;
	cusparseStatus_t stat_t;
	clock_t begin,end;

	stat_t = cusparseCreate(&cusparseHandle);
	stat_t = cusparseCreateMatDescr(&mat_descr);
	stat_t = cusparseSetMatIndexBase(mat_descr,CUSPARSE_INDEX_BASE_ZERO);
	stat_t = cusparseSetMatType(mat_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error!\n");
	}

	int *dev_row_offset,*dev_col;
	double *dev_val,*dev_x,*dev_y;
	const double alpha = 1.0;
	const double beta = 0.0;

	dev_row_offset = (int*)malloc_gpu((nrow+1)*sizeof(int));
	dev_col = (int*)malloc_gpu(nnz*sizeof(int));
	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_row_offset,row,(nrow+1)*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_col,col,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	stat_t = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_TRANSPOSE,nrow,ncol,nnz,&alpha,mat_descr,
							dev_val,dev_row_offset,dev_col,dev_x,&beta,dev_y);
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end-begin)/CLOCKS_PER_SEC;

	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error!\n");
	}

	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row_offset);
	free_gpu(dev_col);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);

	stat_t = cusparseDestroyMatDescr(mat_descr);
	stat_t = cusparseDestroy(cusparseHandle);
}

extern "C" void rmbv_csr(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y){
	cusparseHandle_t cusparseHandle = 0;
	cusparseMatDescr_t mat_descr = 0;
	cusparseStatus_t stat_t;
	clock_t begin,end;

	stat_t = cusparseCreate(&cusparseHandle);
	stat_t = cusparseCreateMatDescr(&mat_descr);
	stat_t = cusparseSetMatIndexBase(mat_descr,CUSPARSE_INDEX_BASE_ZERO);
	stat_t = cusparseSetMatType(mat_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error!\n");
	}

	int *dev_row_offset,*dev_col;
	double *dev_val,*dev_x,*dev_y;
	const double alpha = 1.0;
	const double beta = 0.0;

	dev_row_offset = (int*)malloc_gpu((nrow+1)*sizeof(int));
	dev_col = (int*)malloc_gpu(nnz*sizeof(int));
	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
	dev_x = (double*)malloc_gpu(ncol*sizeof(double));
	dev_y = (double*)malloc_gpu(nrow*sizeof(double));

	cudaMemcpy(dev_row_offset,row,(nrow+1)*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_col,col,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	stat_t = cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,nrow,ncol,nnz,&alpha,mat_descr,
							dev_val,dev_row_offset,dev_col,dev_x,&beta,dev_y);
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;

	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error!\n");
	}

	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row_offset);
	free_gpu(dev_col);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);

	stat_t = cusparseDestroyMatDescr(mat_descr);
	stat_t = cusparseDestroy(cusparseHandle);
}

// extern "C" void rmbv_csr(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
// 	int *dev_row,*dev_PB,*dev_PE;
// 	double *dev_val,*dev_x,*dev_y;
// 	clock_t begin,end;
// 	dim3 BlockSize(GPU_BLOCK_SIZE);
// 	dim3 GridSize((nrow+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

// 	dev_row = (int*)malloc_gpu(nnz*sizeof(int));
// 	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
// 	dev_PB = (int*)malloc_gpu(nrow*sizeof(int));
// 	dev_PE = (int*)malloc_gpu(nrow*sizeof(int));
// 	dev_x = (double*)malloc_gpu(ncol*sizeof(double));
// 	dev_y = (double*)malloc_gpu(nrow*sizeof(double));

// 	cudaMemcpy(dev_row,row,nnz*sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
// 	cudaMemcpy(dev_PB,PB,nrow*sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMemcpy(dev_PE,PE,nrow*sizeof(int),cudaMemcpyHostToDevice);
// 	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
// 	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);

// 	begin = clock();
// 	mbv_csr<<<GridSize,BlockSize>>>(dev_row,dev_val,dev_PB,dev_PE,nrow,dev_x,dev_y);
// 	cudaThreadSynchronize();
// 	end = clock();
// 	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
// 	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

// 	free_gpu(dev_row);
// 	free_gpu(dev_val);
// 	free_gpu(dev_PB);
// 	free_gpu(dev_PE);
// 	free_gpu(dev_x);
// 	free_gpu(dev_y);
// }

extern "C" void lmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int *dev_row,*dev_PB,*dev_PE;
	double *dev_val,*dev_x,*dev_y;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((ncol+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_row = (int*)malloc_gpu(nnz*sizeof(int));
	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
	dev_PB = (int*)malloc_gpu(ncol*sizeof(int));
	dev_PE = (int*)malloc_gpu(ncol*sizeof(int));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_row,row,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PB,PB,ncol*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PE,PE,ncol*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	mbv_csc<<<GridSize,BlockSize>>>(dev_row,dev_val,dev_PB,dev_PE,ncol,dev_x,dev_y,'l');
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row);
	free_gpu(dev_val);
	free_gpu(dev_PB);
	free_gpu(dev_PE);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void rmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int *dev_row,*dev_PB,*dev_PE;
	double *dev_val,*dev_x,*dev_y;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((ncol+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_row = (int*)malloc_gpu(nnz*sizeof(int));
	dev_val = (double*)malloc_gpu(nnz*sizeof(double));
	dev_PB = (int*)malloc_gpu(ncol*sizeof(int));
	dev_PE = (int*)malloc_gpu(ncol*sizeof(int));
	dev_x = (double*)malloc_gpu(ncol*sizeof(double));
	dev_y = (double*)malloc_gpu(nrow*sizeof(double));

	cudaMemcpy(dev_row,row,nnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nnz*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PB,PB,ncol*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PE,PE,ncol*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	mbv_csc<<<GridSize,BlockSize>>>(dev_row,dev_val,dev_PB,dev_PE,ncol,dev_x,dev_y,'r');
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row);
	free_gpu(dev_val);
	free_gpu(dev_PB);
	free_gpu(dev_PE);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void lmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y){
	int *dev_diag;
	double *dev_val,*dev_x,*dev_y;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((ncol+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_diag = (int*)malloc_gpu(ndiag*sizeof(int));
	dev_val = (double*)malloc_gpu(ndiag*lda*sizeof(double));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_diag,diag,ndiag*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,ndiag*lda*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);
	
	begin = clock();
	mbv_dia<<<GridSize,BlockSize>>>(dev_diag,dev_val,nrow,ncol,ndiag,lda,dev_x,dev_y,'l');
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_diag);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void rmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y){
	int *dev_diag;
	double *dev_val,*dev_x,*dev_y;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((nrow+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_diag = (int*)malloc_gpu(ndiag*sizeof(int));
	dev_val = (double*)malloc_gpu(ndiag*lda*sizeof(double));
	dev_x = (double*)malloc_gpu(ncol*sizeof(double));
	dev_y = (double*)malloc_gpu(nrow*sizeof(double));

	cudaMemcpy(dev_diag,diag,ndiag*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,ndiag*lda*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);


	mbv_dia<<<GridSize,BlockSize>>>(dev_diag,dev_val,nrow,ncol,ndiag,lda,dev_x,dev_y,'r');
	
	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_diag);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void lmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int bnnz,int lda,double* x,double* y){
	int *dev_row,*dev_PB,*dev_PE;
	double *dev_val,*dev_x,*dev_y;
	int dummy = lda*lda;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((nb+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_row = (int*)malloc_gpu(bnnz*sizeof(int));
	dev_val = (double*)malloc_gpu(dummy*bnnz*sizeof(double));
	dev_PB = (int*)malloc_gpu(nb*sizeof(int));
	dev_PE = (int*)malloc_gpu(nb*sizeof(int));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_row,row,bnnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,bnnz*dummy*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PB,PB,nb*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PE,PE,nb*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	mbv_bsc<<<GridSize,BlockSize>>>(dev_row,dev_val,dev_PB,dev_PE,nrow,ncol,nb,lda,dev_x,dev_y);
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row);
	free_gpu(dev_val);
	free_gpu(dev_PB);
	free_gpu(dev_PE);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void rmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int bnnz,int lda,double* x,double* y){
	int *dev_row,*dev_PB,*dev_PE;
	double *dev_val,*dev_x,*dev_y;
	int dummy = lda*lda;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((nb+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_row = (int*)malloc_gpu(bnnz*sizeof(int));
	dev_val = (double*)malloc_gpu(dummy*bnnz*sizeof(double));
	dev_PB = (int*)malloc_gpu(nb*sizeof(int));
	dev_PE = (int*)malloc_gpu(nb*sizeof(int));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_row,row,bnnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,bnnz*dummy*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PB,PB,nb*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_PE,PE,nb*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();
	mbv_bsc<<<GridSize,BlockSize>>>(dev_row,dev_val,dev_PB,dev_PE,nrow,ncol,nb,lda,dev_x,dev_y);
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row);
	free_gpu(dev_val);
	free_gpu(dev_PB);
	free_gpu(dev_PE);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void rmbv_bsr(int* row,int* col,double* val,int mb,int nb,int nrow,int ncol,int bnnz,int lda,double* x,double* y){
	cusparseHandle_t cusparseHandle = 0;
	cusparseMatDescr_t mat_descr = 0;
	cusparseStatus_t stat_t;
	clock_t begin,end;

	stat_t = cusparseCreate(&cusparseHandle);
	stat_t = cusparseCreateMatDescr(&mat_descr);
	stat_t = cusparseSetMatIndexBase(mat_descr,CUSPARSE_INDEX_BASE_ZERO);
	stat_t = cusparseSetMatType(mat_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error1!\n");
	}

	int *dev_row_offset,*dev_col;
	double *dev_val,*dev_x,*dev_y;
	const double alpha = 1.0;
	const double beta = 0.0;
	int dummy = lda*lda;

	dev_row_offset = (int*)malloc_gpu((mb+1)*sizeof(int));
	dev_col = (int*)malloc_gpu(bnnz*sizeof(int));
	dev_val = (double*)malloc_gpu(bnnz*dummy*sizeof(double));
	dev_x = (double*)malloc_gpu(nb*lda*sizeof(double));
	dev_y = (double*)malloc_gpu(mb*lda*sizeof(double));

	cudaMemcpy(dev_row_offset,row,(mb+1)*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_col,col,bnnz*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,bnnz*dummy*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);

	begin = clock();	
	stat_t = cusparseDbsrmv(cusparseHandle,CUSPARSE_DIRECTION_ROW,CUSPARSE_OPERATION_NON_TRANSPOSE,mb,nb,bnnz,&alpha,mat_descr,
							dev_val,dev_row_offset,dev_col,lda,dev_x,&beta,dev_y);
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;

	if(CUSPARSE_STATUS_SUCCESS != stat_t){
		printf("Error!\n");
	}

	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_row_offset);
	free_gpu(dev_col);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);

	stat_t = cusparseDestroyMatDescr(mat_descr);
	stat_t = cusparseDestroy(cusparseHandle);
}

extern "C" void lmbv_bdi(int*bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y){
	int *dev_bdiag;
	double *dev_val,*dev_x,*dev_y;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((nb+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_bdiag = (int*)malloc_gpu(nbdiag*sizeof(int));
	dev_val = (double*)malloc_gpu(nbdiag*blda*mm*nn*sizeof(double));
	dev_x = (double*)malloc_gpu(nrow*sizeof(double));
	dev_y = (double*)malloc_gpu(ncol*sizeof(double));

	cudaMemcpy(dev_bdiag,bdiag,nbdiag*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nbdiag*blda*mm*nn*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,nrow*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,ncol*sizeof(double),cudaMemcpyHostToDevice);
	
	begin = clock();
	mbv_bdi<<<GridSize,BlockSize>>>(dev_bdiag,dev_val,nrow,ncol,nbdiag,blda,mb,nb,mm,nn,dev_x,dev_y,'l');
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,ncol*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_bdiag);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);
}

extern "C" void rmbv_bdi(int*bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y){
	int *dev_bdiag;
	double *dev_val,*dev_x,*dev_y;
	clock_t begin,end;
	dim3 BlockSize(GPU_BLOCK_SIZE);
	dim3 GridSize((mb+GPU_BLOCK_SIZE-1)/GPU_BLOCK_SIZE);

	dev_bdiag = (int*)malloc_gpu(nbdiag*sizeof(int));
	dev_val = (double*)malloc_gpu(nbdiag*blda*mm*nn*sizeof(double));
	dev_x = (double*)malloc_gpu(ncol*sizeof(double));
	dev_y = (double*)malloc_gpu(nrow*sizeof(double));

	cudaMemcpy(dev_bdiag,bdiag,nbdiag*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_val,val,nbdiag*blda*mm*nn*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_x,x,ncol*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y,nrow*sizeof(double),cudaMemcpyHostToDevice);
	
	begin = clock();
	mbv_bdi<<<GridSize,BlockSize>>>(dev_bdiag,dev_val,nrow,ncol,nbdiag,blda,mb,nb,mm,nn,dev_x,dev_y,'r');
	cudaThreadSynchronize();
	end = clock();
	compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
	cudaMemcpy(y,dev_y,nrow*sizeof(double),cudaMemcpyDeviceToHost);

	free_gpu(dev_bdiag);
	free_gpu(dev_val);
	free_gpu(dev_x);
	free_gpu(dev_y);
}