#include "sp_io.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>


dsp_linknode* read_coo_matrix(const char *mm_filename)
{
    FILE * fid;
    dsp_linknode* dsp_A;
    MM_typecode matcode;
    int num_rows, num_cols, num_nonzeros;
    int i;

    fid = fopen(mm_filename, "r");
    if (fid == NULL){
        printf("Unable to open file %s\n", mm_filename);
        exit(EXIT_FAILURE);
    }

    if (mm_read_banner(fid, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(EXIT_FAILURE);
    }

    if (!mm_is_valid(matcode)){
        printf("Invalid Matrix Market file.\n");
        exit(EXIT_FAILURE);
    }

    if (!((mm_is_real(matcode) || mm_is_integer(matcode) || mm_is_pattern(matcode)) && mm_is_coordinate(matcode) && mm_is_sparse(matcode) ) ){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        printf("Only sparse real-valued or pattern coordinate matrices are supported\n");
        exit(EXIT_FAILURE);
    }

    if ( mm_read_mtx_crd_size(fid,&num_rows,&num_cols,&num_nonzeros) !=0){
        exit(EXIT_FAILURE);
    }

    int* I = (int *) aligned_malloc(num_nonzeros * sizeof(int));
    int* J = (int *) aligned_malloc(num_nonzeros * sizeof(int));
    double* V = (double*) aligned_malloc(num_nonzeros * sizeof(double));
    memset(I,0,num_nonzeros*sizeof(int));
    memset(J,0,num_nonzeros*sizeof(int));
    memset(V,0.0,num_nonzeros*sizeof(double));

    printf("Reading sparse matrix from file (%s):",mm_filename);
    fflush(stdout);

    if (mm_is_pattern(matcode)){
        // pattern matrix defines sparsity pattern, but not values
        for( i = 0; i < num_nonzeros; i++ ){
            assert(fscanf(fid, " %d %d \n", &I[i], &J[i]) == 2);
            I[i]--;      //adjust from 1-based to 0-based indexing
            J[i]--;
            V[i] = 1.0;  //use value 1.0 for all nonzero entries 
        }
    } else if (mm_is_real(matcode) || mm_is_integer(matcode)){
        for( i = 0; i < num_nonzeros; i++ ){
            int row,col;
            double val;  // always read in a double and convert later if necessary
            
            assert(fscanf(fid, " %d %d %lf \n", &row, &col, &val) == 3);

            I[i] = (int)row - 1; 
            J[i] = (int)col - 1;
            V[i] = (double)val;
        }
    } else {
        printf("Unrecognized data type\n");
        exit(1);
    }

    fclose(fid);
    printf(" done\n");

    if( mm_is_symmetric(matcode) ){ //duplicate off diagonal entries
        int off_diagonals = 0;
        for( i = 0; i < num_nonzeros; i++ ){
            if( I[i] != J[i] )
                off_diagonals++;
        }

        int true_nonzeros = 2*off_diagonals + (num_nonzeros - off_diagonals);

        int* new_I = (int*)aligned_malloc(true_nonzeros*sizeof(int));
        int* new_J = (int*)aligned_malloc(true_nonzeros*sizeof(int));
        double* new_V = (double*)aligned_malloc(true_nonzeros*sizeof(double));
        memset(new_I,0,true_nonzeros*sizeof(int));
        memset(new_J,0,true_nonzeros*sizeof(int));
        memset(new_V,0.0,true_nonzeros*sizeof(double));

        int ptr = 0;
        for( i = 0; i < num_nonzeros; i++ ){
            if( I[i] != J[i] ){
                new_I[ptr] = I[i];  new_J[ptr] = J[i];  new_V[ptr] = V[i];
                ptr++;
                new_J[ptr] = I[i];  new_I[ptr] = J[i];  new_V[ptr] = V[i];
                ptr++;
            } else {
                new_I[ptr] = I[i];  new_J[ptr] = J[i];  new_V[ptr] = V[i];
                ptr++;
            }
        }       
        aligned_free(I); aligned_free(J); aligned_free(V);
        I = new_I; J = new_J; V = new_V;
        num_nonzeros = true_nonzeros;
    } //end symmetric case
    int result = -1;
    dsp_A = duscr_coo(num_rows,num_cols,V,num_nonzeros,I,J,num_nonzeros,1088,&result);
    if(0 == result){
      aligned_free(I);
      aligned_free(J);
      aligned_free(V);
      return dsp_A;
    }else{
      aligned_free(I);
      aligned_free(J);
      aligned_free(V);
      return NULL;
    }
}