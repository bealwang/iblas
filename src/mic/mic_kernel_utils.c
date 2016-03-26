#include "mic/mic_kernel_utils.h"
#include <stdio.h>
#include <time.h>

extern double compute_time;

__attribute__((target(mic))) int mic_dtn(int n,int min_n){
  int max_tn = n/min_n;
  int tn = max_tn>BEST_MIC_THREADS?BEST_MIC_THREADS:max_tn;
  if(tn < 1){
    tn = 1;
  }
  return tn;
}

void lmbv_coo(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y){
	int i;
	#pragma offload target(mic) \
	in(row:length(nnz)) \
	in(col:length(nnz)) \
	in(val:length(nnz)) \
	in(x:length(nrow)) \
	inout(y:length(ncol))
	for(i=0;i<nnz;i++){
		y[col[i]] += val[i]*x[row[i]];
	}
}

void rmbv_coo(int* row,int* col,double* val,int nrow,int ncol,int nnz,double* x,double* y){
	int i;
	#pragma offload target(mic) \
	in(row:length(nnz)) \
	in(col:length(nnz)) \
	in(val:length(nnz)) \
	in(x:length(ncol)) \
	inout(y:length(nrow))
	for(i=0;i<nnz;i++){
		y[row[i]] += val[i]*x[col[i]];
	}
}

void lmbv_csr(int* col,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int i;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(col:length(nnz) align(64) alloc_if(1) free_if(0)) \
	in(val:length(nnz) align(64) alloc_if(1) free_if(0)) \
	in(PB:length(nrow) align(64) alloc_if(1) free_if(0)) \
	in(PE:length(nrow) align(64) alloc_if(1) free_if(0)) \
	in(x:length(nrow) align(64) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) align(64) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(i=0;i<ncol;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(col) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for (i=0;i<nrow;i++){
        int row_begin = PB[i];
        int row_end = PE[i];
        int j;
        for (j=row_begin;j<row_end;j++){
          y[col[j]]+=val[j]*x[i];
        }
    }
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(col:length(nnz) align(64) alloc_if(0) free_if(1)) \
	nocopy(val:length(nnz) align(64) alloc_if(0) free_if(1)) \
	nocopy(PB:length(nrow) align(64) alloc_if(0) free_if(1)) \
	nocopy(PE:length(nrow) align(64) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) align(64) alloc_if(0) free_if(1)) \
	out(y:length(ncol) align(64) alloc_if(0) free_if(1))
	{}
}

void rmbv_csr(int* col,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int i;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(col:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(nnz) alloc_if(1) free_if(0)) \
	in(PB:length(nrow) alloc_if(1) free_if(0)) \
	in(PE:length(nrow) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(i=0;i<nrow;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(col) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for (i=0;i<nrow;i++){
        int row_begin = PB[i];
        int row_end = PE[i];
        int j;
        for (j=row_begin;j<row_end;j++){
          y[i]+=val[j]*x[col[j]];
        }
    }
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(col:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(PB:length(nrow) alloc_if(0) free_if(1)) \
	nocopy(PE:length(nrow) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int j;
	
	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(nnz) alloc_if(1) free_if(0)) \
	in(PB:length(ncol) alloc_if(1) free_if(0)) \
	in(PE:length(ncol) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(j=0;j<ncol;j++){
		y[j] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for (j=0;j<ncol;j++){
        int col_begin = PB[j];
        int col_end = PE[j];
        int i;
        double sum = 0.0;
        for (i=col_begin;i<col_end;i++){
          sum += val[i]*x[row[i]];
        }
        y[j] = sum;
    }
	
    #pragma offload target(mic) \
	nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(PB:length(ncol) alloc_if(0) free_if(1)) \
	nocopy(PE:length(ncol) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_csc(int* row,double* val,int* PB,int* PE,int nrow,int ncol,int nnz,double* x,double* y){
	int j;

	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(nnz) alloc_if(1) free_if(0)) \
	in(PB:length(ncol) alloc_if(1) free_if(0)) \
	in(PE:length(ncol) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(j=0;j<nrow;j++){
		y[j] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for (j=0;j<ncol;j++){
        int col_begin = PB[j];
        int col_end = PE[j];
        int i;
        for (i=col_begin;i<col_end;i++){
          y[row[i]]+=val[i]*x[j];
        }
    }
	
    #pragma offload target(mic) \
	nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(PB:length(ncol) alloc_if(0) free_if(1)) \
	nocopy(PE:length(ncol) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y){
	int i;
	int val_dia = ndiag*lda;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(diag:length(ndiag) alloc_if(1) free_if(0)) \
	in(val:length(val_dia) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(i=0;i<ncol;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(diag) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for (i=0; i<ncol; ++i) {
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

			if ((i >= jstart) && (i < jend)){
				sum += val[i+v_offset] * x[i+I_offset];
			}
			if(i >= jend){
				continue;
			}
		}
		y[i] = sum;
	}
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(diag:length(ndiag) alloc_if(0) free_if(1)) \
	nocopy(val:length(val_dia) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_dia(int* diag,double* val,int nrow,int ncol,int ndiag,int lda,double* x,double* y){
	int i;
	int val_dia = ndiag*lda;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(diag:length(ndiag) alloc_if(1) free_if(0)) \
	in(val:length(val_dia) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(i=0;i<nrow;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(diag) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for (i=0; i<nrow; ++i) {
		double sum = 0.0;
		int j;
		for (j=0; j<ndiag; ++j){
			int offset = diag[j];
			int istart = offset >= 0?0:(-offset);
			int jstart = offset >= 0?offset:0;
			int N = 0; int v_offset = 0; int J_offset = 0;

			if((nrow-istart) > (ncol-jstart)){
			    N = ncol - jstart;
			    J_offset = jstart - istart;
			    v_offset = j*lda - istart;
			}else{
				N = nrow - istart;
				J_offset = jstart - istart;
				v_offset = j*lda + lda - nrow;
			}

			int iend = istart + N;

			if ((i >= istart) && (i < iend)){
				sum += val[i+v_offset] * x[i+J_offset];
			}
			if(i >= iend){
				continue;
			}
		}
		y[i] = sum;
	}
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(diag:length(ndiag) alloc_if(0) free_if(1)) \
	nocopy(val:length(val_dia) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_bco(int* row,int* col,double* val,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int i;
	int dummy = lda*lda*nnz;

	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(col:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(i=0;i<ncol;i++){
		y[i] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(col) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nnz,MIN_ITERATOR_NUM))
	for(i=0;i<nnz;i++){
        int p;
        for (p=0;p<lda;p++){
	        int q;
	        for (q=0;q<lda;q++)
	        {
	            if((col[i]*lda+q < ncol) && (row[i]*lda+p < nrow)){
	            	y[col[i]*lda+q]+=val[i*lda*lda+p*lda+q]*x[row[i]*lda+p];
	            }
	        }
        }
    }
	
    #pragma offload target(mic) \
    nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(col:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_bco(int* row,int* col,double* val,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int i;
	int dummy = lda*lda*nnz;

	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(col:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(i=0;i<nrow;i++){
		y[i] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(col) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nnz,MIN_ITERATOR_NUM))
	for(i=0;i<nnz;i++){
        int p;
        for (p=0;p<lda;p++){
	        int q;
	        for(q=0;q<lda;q++){
	          if((row[i]*lda+p < nrow) && (col[i]*lda+q < ncol)){
	            y[row[i]*lda+p]+=val[i*lda*lda+p*lda+q]*x[col[i]*lda+q];
	        }
        }
      }
    }
	
    #pragma offload target(mic) \
    nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(col:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_bsr(int* col,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int i;
	int dummy = lda*lda*nnz;

	#pragma offload target(mic) \
	in(col:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(PB:length(mb) alloc_if(1) free_if(0)) \
	in(PE:length(mb) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(i=0;i<ncol;i++){
		y[i] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(col) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(mb,MIN_ITERATOR_NUM))
	for(i=0;i<mb;i++){
        int row_begin = PB[i];
        int row_end = PE[i];
        int j;
        for(j=row_begin;j<row_end;j++)
        {
            int p;
            for (p=0;p<lda;p++){
	            int q;
	            for (q=0;q<lda;q++){
	                if((col[j]*lda+q < ncol) && (i*lda+p < nrow)){
	                    y[col[j]*lda+q]+=val[j*lda*lda+p*lda+q]*x[i*lda+p];
	                }
	            }
            }
        }
    }
	
    #pragma offload target(mic) \
	nocopy(col:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(PB:length(mb) alloc_if(0) free_if(1)) \
	nocopy(PE:length(mb) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_bsr(int* col,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int i;
	int dummy = lda*lda*nnz;

	#pragma offload target(mic) \
	in(col:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(PB:length(mb) alloc_if(1) free_if(0)) \
	in(PE:length(mb) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(i=0;i<nrow;i++){
		y[i] = 0.0;
	}

	#pragma offload target(mic) \
	nocopy(col) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(mb,MIN_ITERATOR_NUM))
	for(i=0;i<mb;i++){
        int row_begin = PB[i];
        int row_end = PE[i];
        int j;
        for(j=row_begin;j<row_end;j++)
        {
            int p;
            for (p=0;p<lda;p++){
	            int q;
	            for (q=0;q<lda;q++){
	                if((col[j]*lda+q < ncol) && (i*lda+p < nrow)){
	                    y[i*lda+p]+=val[j*lda*lda+p*lda+q]*x[col[j]*lda+q];
	                }
	            }
            }
        }
    }
	
    #pragma offload target(mic) \
	nocopy(col:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(PB:length(mb) alloc_if(0) free_if(1)) \
	nocopy(PE:length(mb) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int j;
	int dummy = lda*lda*nnz;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(PB:length(nb) alloc_if(1) free_if(0)) \
	in(PE:length(nb) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(j=0;j<ncol;j++){
		y[j] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nb,MIN_ITERATOR_NUM))
	for(j=0;j<nb;j++){
        int col_begin = PB[j];
        int col_end = PE[j];
        int i;
        for(i=col_begin;i<col_end;i++)
        {
            int p;
            for (p=0;p<lda;p++){
	            int q;
	            for (q=0;q<lda;q++){
	                if((j*lda+q < ncol) && (row[i]*lda+p < nrow)){
	                    y[j*lda+q] += val[i*lda*lda+p*lda+q]*x[row[i]*lda+p];
	                }
	            }
            }
        }
    }
    end = clock();
	compute_time += (double)(end-begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(PB:length(nb) alloc_if(0) free_if(1)) \
	nocopy(PE:length(nb) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_bsc(int* row,double* val,int* PB,int* PE,int mb,int nb,int nrow,int ncol,int nnz,int lda,double* x,double* y){
	int j;
	int dummy = lda*lda*nnz;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(row:length(nnz) alloc_if(1) free_if(0)) \
	in(val:length(dummy) alloc_if(1) free_if(0)) \
	in(PB:length(nb) alloc_if(1) free_if(0)) \
	in(PE:length(nb) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(j=0;j<nrow;j++){
		y[j] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(row) \
	nocopy(val) \
	nocopy(PB) \
	nocopy(PE) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nb,MIN_ITERATOR_NUM))
	for(j=0;j<nb;j++){
        int col_begin = PB[j];
        int col_end = PE[j];
        int i;
        for(i=col_begin;i<col_end;i++)
        {
            int p;
            for (p=0;p<lda;p++){
	            int q;
	            for (q=0;q<lda;q++){
	                if((j*lda+q < ncol) && (row[i]*lda+p < nrow)){
	                    y[row[i]*lda+p] += val[i*lda*lda+p*lda+q]*x[j*lda+q];
	                }
	            }
            }
        }
    }
    end = clock();
	compute_time += (double)(end-begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(row:length(nnz) alloc_if(0) free_if(1)) \
	nocopy(val:length(dummy) alloc_if(0) free_if(1)) \
	nocopy(PB:length(nb) alloc_if(0) free_if(1)) \
	nocopy(PE:length(nb) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}

void lmbv_bdi(int* bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y){
	int i;
	int dummy = mm*nn;
	int val_bdi = nbdiag*blda*dummy;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(bdiag:length(nbdiag) alloc_if(1) free_if(0)) \
	in(val:length(val_bdi) alloc_if(1) free_if(0)) \
	in(x:length(nrow) alloc_if(1) free_if(0)) \
	nocopy(y:length(ncol) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(ncol,MIN_ITERATOR_NUM))
	for(i=0;i<ncol;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(bdiag) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(nb,MIN_ITERATOR_NUM))
	for (i=0; i<nb; ++i) {
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
	  
	  if ((i >= jstart) && (i < jend)){
	    int p;
	    for (p=0;p<nn;p++){
	      int q;
	      for (q=0;q<nn;q++){
	        if((i*nn+q < ncol) && ((i+I_offset)*mm+p < nrow)){
	          y[i*nn+q] += val[(i+v_offset)*mm*nn+nn*p+q] * x[(i+I_offset)*mm+p];
	        }
	      }
	    }
	  }
	  if(i >= jend){
	    continue;
	  }
	}
	}
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(bdiag:length(nbdiag) alloc_if(0) free_if(1)) \
	nocopy(val:length(val_bdi) alloc_if(0) free_if(1)) \
	nocopy(x:length(nrow) alloc_if(0) free_if(1)) \
	out(y:length(ncol) alloc_if(0) free_if(1))
	{}
}

void rmbv_bdi(int* bdiag,double* val,int nrow,int ncol,int nbdiag,int blda,int mb,int nb,int mm,int nn,double* x,double* y){
	int i;
	int dummy = mm*nn;
	int val_bdi = nbdiag*blda*dummy;
	clock_t begin,end;

	#pragma offload target(mic) \
	in(bdiag:length(nbdiag) alloc_if(1) free_if(0)) \
	in(val:length(val_bdi) alloc_if(1) free_if(0)) \
	in(x:length(ncol) alloc_if(1) free_if(0)) \
	nocopy(y:length(nrow) alloc_if(1) free_if(0))
	#pragma omp parallel for num_threads(mic_dtn(nrow,MIN_ITERATOR_NUM))
	for(i=0;i<nrow;i++){
		y[i] = 0.0;
	}

	begin = clock();
	#pragma offload target(mic) \
	nocopy(bdiag) \
	nocopy(val) \
	nocopy(x) \
	nocopy(y)
	#pragma omp parallel for num_threads(mic_dtn(mb,MIN_ITERATOR_NUM))
	for (i=0; i<mb; ++i) {
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
	  
	  if ((i >= istart) && (i < iend)){
	    int p;
	    for (p=0;p<mm;p++){
	      int q;
	      for (q=0;q<nn;q++){
	        if((i*mm+p <nrow) && ((i+J_offset)*nn+q < ncol)){
	          y[i*mm+p] += val[(i+v_offset)*mm*nn+nn*p+q] * x[(i+J_offset)*nn+q];
	        }
	      }
	    }
	  }
	  if(i >= iend){
	    continue;
	  }
	}
	}
    end = clock();
    compute_time += (double)(end - begin)/CLOCKS_PER_SEC;
	
    #pragma offload target(mic) \
	nocopy(bdiag:length(nbdiag) alloc_if(0) free_if(1)) \
	nocopy(val:length(val_bdi) alloc_if(0) free_if(1)) \
	nocopy(x:length(ncol) alloc_if(0) free_if(1)) \
	out(y:length(nrow) alloc_if(0) free_if(1))
	{}
}