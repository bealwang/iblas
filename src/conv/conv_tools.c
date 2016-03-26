#include "conv/conv_tools.h"
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

void dpre_usconv_coo2csr(int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int* PB,int* PE)
{
  int i,cumsum,last;

  int* Aj = (int*)aligned_malloc(sizeof(int)*num_nonzeros);
  int* temp_ind = (int*)aligned_malloc(sizeof(int)*num_rows);
  double* Ax = (double*)aligned_malloc(sizeof(double)*num_nonzeros);
  memset(Aj,0,sizeof(int)*num_nonzeros);
  memset(temp_ind,0,sizeof(int)*num_rows);
  memset(Ax,0.0,sizeof(double)*num_nonzeros);

  for (i = 0; i < num_nonzeros; i++) 
    PB[(*rows)[i]]++;

  //cumsum the nnz per row to get Bp[]
  for(i = 0,cumsum = 0;i < num_rows; i++){     
      int temp = PB[i];
      PB[i] = cumsum;
      cumsum += temp;
  }

  for(i = 0;i < num_rows-1; i++){
      PE[i] = PB[i+1];
  }
  PE[i] = num_nonzeros;

  memcpy(temp_ind,PB,num_rows*sizeof(int));
  //write Aj,Ax into Bj,Bx
  for(i = 0;i < num_nonzeros; i++){
      int row  = (*rows)[i];
      int dest = temp_ind[row];
      Aj[dest] = cols[i];
      Ax[dest] = (*data)[i];
      temp_ind[row]++;
  }

  aligned_free(temp_ind);
  aligned_free(*rows);
  aligned_free(*data);

  *rows = Aj;
  *data = Ax;
}

void dpre_usconv_coo2csc(int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int* PB,int* PE)
{
  int i,cumsum,last;

  int* Aj = (int*)aligned_malloc(sizeof(int)*num_nonzeros);
  int* temp_ind = (int*)aligned_malloc(sizeof(int)*num_cols);
  double* Ax = (double*)aligned_malloc(sizeof(double)*num_nonzeros);
  memset(Aj,0,sizeof(int)*num_nonzeros);
  memset(temp_ind,0,sizeof(int)*num_cols);
  memset(Ax,0.0,sizeof(double)*num_nonzeros);

  for (i = 0; i < num_nonzeros; i++)
    PB[cols[i]]++;

  //cumsum the nnz per row to get Bp[]
  for(i = 0,cumsum = 0;i < num_cols; i++){
      int temp = PB[i];
      PB[i] = cumsum;
      cumsum += temp;
  }

  for(i = 0;i < num_cols-1; i++){
      PE[i] = PB[i+1];
  }
  PE[i] = num_nonzeros;

  memcpy(temp_ind,PB,num_cols*sizeof(int));
  //write Aj,Ax into Bj,Bx
  for(i = 0;i < num_nonzeros; i++){
      int col  = cols[i];
      int dest = temp_ind[col];
      Aj[dest] = (*rows)[i];
      Ax[dest] = (*data)[i];
      temp_ind[col]++;
  }

  aligned_free(temp_ind);
  aligned_free(*rows);
  aligned_free(*data);

  *rows = Aj;
  *data = Ax;
}

void dpre_usconv_coo2dia  (int m,int n,double** VAL,int n_VAL,int** INDX,int* JNDX,int* LDA,int* NDIAG)
{
  *NDIAG=0;
  *LDA=m<n?m:n;
  int num_rows,num_cols,num_nonzeros;
  int i,ii,jj,offset,map_index,complete_ndiags;
  size_t VAL_DIA_size = 0;
  int unmarked = -1;
  num_rows = m;
  num_cols = n;
  num_nonzeros = n_VAL;
  complete_ndiags = 0;
  int* diag_map = (int*)aligned_malloc(sizeof(int)*(num_rows+num_cols));     //mark the diags
  int* diag_offset = (int*)aligned_malloc(sizeof(int)*(num_cols+num_rows));
  memset(diag_map,unmarked,sizeof(int)*(num_rows+num_cols));
  memset(diag_offset,unmarked,sizeof(int)*(num_rows+num_cols));

  for(i=0;i<num_nonzeros;i++){
    ii = (*INDX)[i];
    jj = JNDX[i];
    map_index = num_rows-ii+jj;            //be used to find the same diag
    if(diag_map[map_index] == unmarked){        
      diag_map[map_index] = complete_ndiags;
      diag_offset[map_index] = jj - ii;    //get index of diags
      complete_ndiags++;                 //number of diags
    }
  }

  double max_fill = NUM_DIAGS_LIMIT_PARA;
  int dia_max_diags = (max_fill*num_nonzeros)/num_rows + 1;
  if(complete_ndiags > dia_max_diags)
    {                                                                      
       printf("\tNumber of diagonals (%d) excedes limit (%d)\n",complete_ndiags, dia_max_diags);
       aligned_free(diag_map);
       aligned_free(diag_offset);                             
       return;
   }

  *NDIAG = complete_ndiags;
  VAL_DIA_size=(size_t)(*NDIAG)*(*LDA);
  int* IDIAG=(int*)aligned_malloc(sizeof(int)*(complete_ndiags));
  double* VAL_DIA=(double*)aligned_malloc(sizeof(double)*VAL_DIA_size);
  memset(IDIAG,0,sizeof(int)*(complete_ndiags));
  memset(VAL_DIA,0.0,sizeof(double)*VAL_DIA_size);

  for(i=0;i<num_rows + num_cols;i++){
    if(diag_map[i] != unmarked){
      int idiag_ind = diag_map[i];
      IDIAG[idiag_ind] = diag_offset[i];     //offset of diags
    }
  }
  
  for(i=0;i<num_nonzeros;i++){    //get values of every diag
    int ii = (*INDX)[i];
    int jj = JNDX[i];
    int map_index = num_rows-ii+jj;
    int diag = diag_map[map_index];
    if(diag_offset[map_index] >=0){
      offset = ii;
    }else{
      int istart = -diag_offset[map_index];
      if((num_rows-istart) >= num_cols){
        offset = ii - istart;
      }else{
        offset = (*LDA) + ii - num_rows;
      }
      
    }
    VAL_DIA[(size_t)(diag*(*LDA)+offset)] = (*VAL)[i];
  }

  aligned_free(diag_map);
  aligned_free(diag_offset);
  aligned_free(*VAL);
  aligned_free(*INDX);
  //把结果的值拷贝到VAL，INDX中
  //其中VAL_DIA中存放斜线中的值
  //INDX中存放斜线号与斜线的偏移值
  *VAL = VAL_DIA;
  *INDX = IDIAG;
}
void dpre_usconv_coo2bco  (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int row_block_size,int col_block_size,int* bnnz,int* mb,int* kb){
    int num_block_rows = (m+row_block_size-1) / row_block_size;
    int num_block_cols = (n+col_block_size-1) / col_block_size;
    int num_rows_left = m % row_block_size;
    int block_size = row_block_size * col_block_size;

    int *bAi;
    int *bAj;
    double *bAx;
    
    //tmp variables
    int num_blocks = 0;
    int** block_count;
    int i, j,i0, j0, k, I, J, tmp_num_rows_left;
    //block_count is used to mark the block that have been traversed
    block_count = (int**)malloc2d(num_block_rows,num_block_cols,sizeof(int));
    initarray(block_count,num_block_rows,num_block_cols);
    //initarray(block_count,num_block_rows,num_block_cols);
    tmp_num_rows_left = num_rows_left == 0 ? row_block_size : num_rows_left;
    //Phase I: Count the exact number of new blocks to create. 

    //#pragma omp parallel for num_threads(dtn(n_VAL,MIN_ITERATOR_NUM))
    for(i=0;i<n_VAL;i++){
      I = (*INDX)[i] / row_block_size;
      J = (*JNDX)[i] / col_block_size;
      if(block_count[I][J] == 0){
          num_blocks ++;
          block_count[I][J] ++;
      }
    }
    //End of Phase I
    *bnnz = num_blocks;
    *mb = num_block_rows;
    *kb = num_block_cols;
    
    bAi = (int*)aligned_malloc(sizeof(int)*num_blocks);
    bAj = (int*)aligned_malloc(sizeof(int)*num_blocks);
    bAx = (double*)aligned_malloc(sizeof(double)*num_blocks*block_size);
    double *blocks = (double*)aligned_malloc(sizeof(double)*num_block_rows*num_block_cols*block_size);
    memset (blocks, 0, sizeof(double)*num_block_rows*num_block_cols*block_size);
    memset (bAx, 0, sizeof(double)*num_blocks*block_size);
    //Phase II: Copy values of all blocks.    
    int nnzIndex = 0;
    #pragma omp parallel for num_threads(dtn(n_VAL,MIN_ITERATOR_NUM)) private(i,I,i0,j,J,j0)
    for(k=0;k<n_VAL;k++){
        i = (*INDX)[k];
        I = i / row_block_size;
        i0 =i - I*row_block_size;
        j = (*JNDX)[k];
        J = j / col_block_size;
        j0 =j - J * col_block_size;
        blocks[(I*num_block_cols+J)*block_size + i0*col_block_size + j0] = (*VAL)[k];
      }

    //#pragma omp parallel for num_threads(dtn(n_VAL,MIN_ITERATOR_NUM))
    for(k=0;k<n_VAL;k++){
        i = (*INDX)[k];
        I = i / row_block_size;
        j = (*JNDX)[k];
        J = j / col_block_size;
        if(block_count[I][J] > 0){
            memcpy(bAx+nnzIndex*block_size, blocks+(I*num_block_cols+J)*block_size, block_size*sizeof(double));
            bAi[nnzIndex] = I;
            bAj[nnzIndex] = J;
            nnzIndex ++;
            block_count[I][J] = 0;
        }
      }
      //End of Phase II

    aligned_free(block_count);
    aligned_free(blocks);
    aligned_free(*INDX);
    aligned_free(*JNDX);
    aligned_free(*VAL);

    *INDX = bAi;             //copy the result to BCO_format matrix
    *JNDX = bAj;
    *VAL = bAx;
    // printf("\n");
    // for(i=0;i<2;i++){
    //   printf("%d\t%d",bAi[i],bAj[i]);
    // }
    // printf("\n");
    // printf("\n");
    // for(i=0;i<num_blocks*block_size;i++){
    //   printf("%lf\t",bAx[i]);
    // }
    // printf("\n");
}

void dpre_usconv_bco2coo  (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int* nnz,int bnnz,int lb,int mb,int kb){
    int num_block_rows = mb;
    int num_block_cols = kb;
    int row_block_size = lb;
    int block_size = lb*lb;
    int num_blocks = bnnz;
    int i,j,k,ii,jj,num_nonzeros,ind,row_ind,col_ind,val_ind,VAL_ind;
    int *Ai,*Aj;
    double *Ax,val_val;

    num_nonzeros = 0;
    //Phase I : get num of nonzeros
    for(i=0;i<n_VAL;i++){
      if((*VAL)[i] !=0){
        num_nonzeros++;
      }
    }
    //end of Phase I
    *nnz = num_nonzeros;
    Ai = (int*)aligned_malloc(num_nonzeros*sizeof(int));
    Aj = (int*)aligned_malloc(num_nonzeros*sizeof(int));
    Ax = (double*)aligned_malloc(num_nonzeros*sizeof(double));
    memset(Ai,0,num_nonzeros*sizeof(int));
    memset(Aj,0,num_nonzeros*sizeof(int));
    memset(Ax,0.0,num_nonzeros*sizeof(double));

    row_ind = col_ind = val_ind = ind = 0;
    //Phase II :transfer VAL to two-dimensional array
    for(i=0;i<num_blocks;i++){
      ii = (*INDX)[i];
      jj = (*JNDX)[i];
      for(j=0;j<row_block_size;j++){
        for(k=0;k<row_block_size;k++){
          VAL_ind = i*block_size + j*row_block_size +k;
          val_val = (*VAL)[VAL_ind];
          if(val_val != 0){
            row_ind = ii*row_block_size + j;
            col_ind = jj*row_block_size + k;
            val_ind = row_ind*n + col_ind;
            Ai[ind] = row_ind;
            Aj[ind] = col_ind;
            Ax[ind] = val_val;
            ind++;
          }
        }
      }
    }
    //end of Phase II

    aligned_free(*INDX);
    aligned_free(*JNDX);
    aligned_free(*VAL);

    *INDX = Ai;             //copy the result to COO_format matrix
    *JNDX = Aj;
    *VAL = Ax;
}

void dpre_usconv_dia2coo (double** VAL_DIA,int** IDIAG,int n_IDIAG ,int** IA2,int LDA,int NNZ,int m)
{
  double *VAL;
  int *INDX,*JNDX;
  int VAL_size,i,k;
  int VAL_ind,IND_ind,VAL_DIA_ind;
  int nz;//number of padded zero

  VAL_size =NNZ;
  VAL=(double*)aligned_malloc(sizeof(double)*VAL_size);
  INDX=(int*)aligned_malloc(sizeof(int)*VAL_size);
  JNDX=(int*)aligned_malloc(sizeof(int)*VAL_size);

  VAL_ind=0;
  IND_ind=0;
  VAL_DIA_ind=0;
  for(i=0;i<n_IDIAG;i++)
   {
      if((*IDIAG)[i]<0){//low diag
          nz=(m+(*IDIAG)[i])>LDA?0:LDA-(m+(*IDIAG)[i]);
          for(k=0;k<LDA;k++)
         {    if((*VAL_DIA)[VAL_DIA_ind]!=0)//convert only the value is not equal to zero
              {
              VAL[VAL_ind]=(*VAL_DIA)[VAL_DIA_ind];
              INDX[IND_ind]=(-(*IDIAG)[i])+k-nz;
              JNDX[IND_ind]=k-nz;

              IND_ind++;
              VAL_ind++;
          }
              VAL_DIA_ind++;

          }
        }
      else if((*IDIAG)[i]==0){ //main diag
          for(k=0;k<LDA;k++)
           {
              if((*VAL_DIA)[VAL_DIA_ind]!=0){
              VAL[VAL_ind]=(*VAL_DIA)[VAL_DIA_ind];
              INDX[IND_ind]=k;
              JNDX[IND_ind]=k;
              IND_ind++;
              VAL_ind++;
           }
              VAL_DIA_ind++;
        }
        }
     else{//high diag
          for(k=0;k<LDA;k++){
       if((*VAL_DIA)[VAL_DIA_ind]!=0) {
       VAL[VAL_ind]=(*VAL_DIA)[VAL_DIA_ind];
       JNDX[IND_ind]=(*IDIAG)[i]+k;
       INDX[IND_ind]=k;
       IND_ind++;
       VAL_ind++;
         }
       VAL_DIA_ind++;
        }
     }
  }

  aligned_free(*VAL_DIA);
  aligned_free(*IDIAG);

  (*VAL_DIA)=VAL;
  (*IDIAG)=INDX;
  (*IA2)=JNDX;
}

void dpre_usconv_bco2bsr (int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int lb,int* PB,int* PE)
{
  int i,cumsum,last;

  int* Aj = (int*)aligned_malloc(sizeof(int)*num_nonzeros);
  int* temp_ind = (int*)aligned_malloc(sizeof(int)*num_rows);
  double* Ax = (double*)aligned_malloc(sizeof(double)*num_nonzeros*lb*lb);
  memset(Aj,0,sizeof(int)*num_nonzeros);
  memset(temp_ind,0,sizeof(int)*num_rows);
  memset(Ax,0.0,sizeof(double)*num_nonzeros*lb*lb);

  for (i = 0; i < num_nonzeros; i++) 
    PB[(*rows)[i]]++;

  //cumsum the nnz per row to get Bp[]
  for(i = 0,cumsum = 0;i < num_rows; i++){     
      int temp = PB[i];
      PB[i] = cumsum;
      cumsum += temp;
  }

  for(i = 0;i < num_rows-1; i++){
      PE[i] = PB[i+1];
  }
  PE[i] = num_nonzeros;

  memcpy(temp_ind,PB,num_rows*sizeof(int));
  //write Aj,Ax into Bj,Bx
  for(i = 0;i < num_nonzeros; i++){
      int row  = (*rows)[i];
      int dest = temp_ind[row];
      Aj[dest] = cols[i];
      memcpy(Ax+dest*lb*lb,(*data)+i*lb*lb,sizeof(double)*lb*lb);
      temp_ind[row]++;
  }

  aligned_free(temp_ind);
  aligned_free(*rows);
  aligned_free(*data);

  *rows = Aj;
  *data = Ax;

}
void dpre_usconv_bco2bsc (int** rows,int* cols,double** data,int num_rows,int num_cols,int num_nonzeros,int lb,int* PB,int* PE)
{

  int i,cumsum,last;

  int* Aj = (int*)aligned_malloc(sizeof(int)*num_nonzeros);
  int* temp_ind = (int*)aligned_malloc(sizeof(int)*num_cols);
  double* Ax = (double*)aligned_malloc(sizeof(double)*num_nonzeros*lb*lb);
  memset(Aj,0,sizeof(int)*num_nonzeros);
  memset(temp_ind,0,sizeof(int)*num_cols);
  memset(Ax,0.0,sizeof(double)*num_nonzeros*lb*lb);

  for (i = 0; i < num_nonzeros; i++)
    PB[cols[i]]++;

  //cumsum the nnz per row to get Bp[]
  for(i = 0,cumsum = 0;i < num_cols; i++){
      int temp = PB[i];
      PB[i] = cumsum;
      cumsum += temp;
  }

  for(i = 0;i < num_cols-1; i++){
      PE[i] = PB[i+1];
  }
  PE[i] = num_nonzeros;

  memcpy(temp_ind,PB,num_cols*sizeof(int));
  //write Aj,Ax into Bj,Bx
  for(i = 0;i < num_nonzeros; i++){
      int col  = cols[i];
      int dest = temp_ind[col];
      Aj[dest] = (*rows)[i];
      memcpy(Ax+dest*lb*lb,(*data)+i*lb*lb,sizeof(double)*lb*lb);
      temp_ind[col]++;
  }

  aligned_free(temp_ind);
  aligned_free(*rows);
  aligned_free(*data);

  *rows = Aj;
  *data = Ax;
}
void dpre_usconv_bco2bdi (int mb,int kb,int lb,double** VAL,int** BINDX,int n_BINDX,int* BJNDX,int* BLDA,int* BNDIAG)
{
  int i,ii,jj,offset,val_offset;
  const int unmarked = -1;
  int num_rows,num_cols,num_nonzeros,map_index,complete_ndiags,VAL_DIA_size,dummy;

  *BLDA=mb<kb?mb:kb;
  (*BNDIAG)=0;
  dummy = lb*lb;
  num_rows = mb;
  num_cols = kb;
  num_nonzeros = n_BINDX;
  complete_ndiags = 0;

  int* diag_map = (int*)aligned_malloc((num_rows+num_cols)*sizeof(int));
  int* diag_offset = (int*)aligned_malloc((num_rows+num_cols)*sizeof(int));
  memset(diag_map,unmarked,(num_cols+num_rows)*sizeof(int));
  memset(diag_offset,unmarked,(num_cols+num_rows)*sizeof(int));

  //Phase I:count the number of diags
  for(i=0;i<num_nonzeros;i++){
    ii = (*BINDX)[i];
    jj = BJNDX[i];
    map_index = num_rows-ii+jj;
    if(diag_map[map_index] == unmarked){
      diag_map[map_index] = complete_ndiags;
      complete_ndiags++;
      diag_offset[map_index] = jj-ii;
    }
  }
  //end of Phase I

  double max_fill = NUM_DIAGS_LIMIT_PARA;
  int bdi_max_diags = (max_fill*num_nonzeros)/mb + 1;
  if(complete_ndiags > bdi_max_diags)
    {                                                                      
       printf("\tNumber of block_diagonals (%d) excedes limit (%d)\n",complete_ndiags, bdi_max_diags);
       aligned_free(diag_map);
       aligned_free(diag_offset);                             
       return;
   }

  *BNDIAG = complete_ndiags;
  VAL_DIA_size=(*BNDIAG)*(*BLDA)*dummy;
  int* BIDIAG=(int*)aligned_malloc((*BNDIAG)*sizeof(int));
  double* VAL_DIA = (double*)aligned_malloc(VAL_DIA_size*sizeof(double));
  //double* VAL_DIA=(double*)new_array(VAL_DIA_size,sizeof(double));
  // if(BIDIAG != NULL && VAL_DIA != NULL){
  //   printf("error\n");
  // }
  memset(BIDIAG,0,(*BNDIAG)*sizeof(int));
  //printf("2\n");
  memset(VAL_DIA,0,VAL_DIA_size*sizeof(double));
  //printf("3\n");
  //Phase II:get offset of every diag
  for(i=0;i<num_rows + num_cols;i++){
    if(diag_map[i] != unmarked){
      BIDIAG[diag_map[i]] = diag_offset[i];
    }
  }
  //end of Phase II
  //printf("begin\n");
  //Phase III:get values of every diag
  for(i=0;i<num_nonzeros;i++){
    ii = (*BINDX)[i];
    jj = BJNDX[i];
    map_index = num_rows-ii+jj;
    int diag = diag_map[map_index];
    if((*BLDA)-(num_rows-ii)>0){
      offset=(*BLDA)-(num_rows-ii);
    }else{
      offset=0;
    }
    val_offset = diag*(*BLDA)*dummy + offset*dummy;
    memcpy(VAL_DIA+val_offset,(*VAL)+i*dummy,dummy*sizeof(double));
  }
  //end of Phase III
  //printf("1\n");
  aligned_free(diag_map);
  aligned_free(diag_offset);
  aligned_free(*VAL);
  aligned_free(*BINDX);

  *VAL=VAL_DIA;
  *BINDX=BIDIAG;
 
  // printf("\nBIDIAG:");
  // for(i=0;i<complete_ndiags;i++){
  //   printf("%d\t",BIDIAG[i]);
  // }
  // printf("\nVAL_DIA:");
  // for(i=0;i<VAL_DIA_size;i++){
  //   printf("%lf\t",VAL_DIA[i]);
  // }
}

void dpre_usconv_bdi2bco (double** VAL_DIA,int n_VAL_DIA,int** BIDIAG,int n_BIDIAG,int** IA2,int BLDA,int* BNNZ,int lb)
{

  double *VAL;
  int *BINDX,*BJNDX;
  int val_val;
  int VAL_size,i,k,VAL_ind,IND_ind;
  int VAL_DIA_ind,dummy,sub_ind,NB_BLOCKS,dummy2;
  bool sub_finder;
  dummy=lb*lb;
  VAL_size =*BNNZ;
  //*******get the VAL BINDX BJNDX size********
  VAL=(double*)aligned_malloc(dummy*VAL_size*sizeof(double));
  BINDX=(int*)aligned_malloc (VAL_size*sizeof(int));
  BJNDX=(int*)aligned_malloc (VAL_size*sizeof(int));
  memset(BINDX,0,VAL_size*sizeof(int));
  memset(BJNDX,0,VAL_size*sizeof(int));
  memset(VAL,0,dummy*VAL_size*sizeof(double));

  VAL_ind=0;
  IND_ind=0;
  VAL_DIA_ind=0;

  for(i=0;i<n_BIDIAG;i++)
    {
      if((*BIDIAG)[i]<0) {
          for(k=0;k<BLDA;k++)
            {
              sub_finder=FALSE;
              for(sub_ind =0;sub_ind<dummy;sub_ind++)
                {
                  val_val=(*VAL_DIA)[dummy*VAL_DIA_ind+sub_ind];
                  if(val_val!=0)
                    {
                      sub_finder=TRUE;
                      break;
                    }
                }
              if(sub_finder)
                {
                  for(sub_ind=0;sub_ind<dummy;sub_ind++)
                    {
                      dummy2=dummy*VAL_ind+sub_ind;
                      VAL[dummy2]=(*VAL_DIA)[VAL_DIA_ind*dummy+sub_ind];
                    }
                  BINDX[IND_ind]=k-(*BIDIAG)[i];
                  BJNDX[IND_ind]=k;
                  IND_ind++;
                  VAL_ind++;
                }
              VAL_DIA_ind++;
            }
        }
      else if((*BIDIAG)[i]==0)
        {
          for(k=0;k<BLDA;k++)
            {
              sub_finder=FALSE;
              for(sub_ind =0;sub_ind<dummy;sub_ind++)
                {
                  val_val=(*VAL_DIA)[dummy*(VAL_DIA_ind)+sub_ind];
                  if(val_val!=0)
                    {
                      sub_finder=TRUE;
                      break;
                    }

                }
              if(sub_finder)
                {
                  for(sub_ind=0;sub_ind<dummy;sub_ind++)
                    {
                      dummy2=dummy*VAL_ind+sub_ind;
                      VAL[dummy2]=(*VAL_DIA)[VAL_DIA_ind*dummy+sub_ind];
                    }
                  BINDX[IND_ind]=k;
                  BJNDX[IND_ind]=k;
                  IND_ind++;
                  VAL_ind++;
                }
              VAL_DIA_ind++;
            }
        }
      else
        {
          for(k=0;k<BLDA;k++){
            sub_finder=FALSE;
            for(sub_ind =0;sub_ind<dummy;sub_ind++)
              {
                val_val=(*VAL_DIA)[dummy*(VAL_DIA_ind)+sub_ind];
                if(val_val!=0)
                  {
                    sub_finder=TRUE;
                    break;
                  }

              }
            if(sub_finder)
              {
                for(sub_ind=0;sub_ind<dummy;sub_ind++)
                  {
                    dummy2=dummy*VAL_ind+sub_ind;
                    VAL[dummy2]=(*VAL_DIA)[VAL_DIA_ind*dummy+sub_ind];
                  }
                BINDX[IND_ind]=k;
                BJNDX[IND_ind]=k+(*BIDIAG)[i];
                IND_ind++;
                VAL_ind++;
              }
            VAL_DIA_ind++;
          }
        }
    }
      
  aligned_free((*VAL_DIA));
  aligned_free((*BIDIAG));

  (*VAL_DIA)=VAL;
  (*BIDIAG)=BINDX;
  (*IA2)=BJNDX;
}

void dpre_usconv_csr2bco  (int m,int n,double** VAL,int n_VAL,int** INDX,int** JNDX,int* PB,int* PE,int row_block_size,int col_block_size,int* bnnz,int* mb,int* kb)
{  
    int num_block_rows = (m+row_block_size-1) / row_block_size;
    int num_block_cols = (n+col_block_size-1) / col_block_size;
    int num_rows_left = m % row_block_size;
    int block_size = row_block_size * col_block_size;

    int *bAi;
    int *bAj;
    double *bAx;
    int num_blocks = 0;
    int *block_count;
    int i, j, j0, k, I, J, di, tmp_num_rows_left;

    block_count = (int*)aligned_malloc(num_block_cols*sizeof(int));
    memset (block_count, 0, sizeof(int)*num_block_cols);
    tmp_num_rows_left = num_rows_left == 0 ? row_block_size : num_rows_left;
      
    //Phase I: Count the exact number of new blocks to create.    
      for(I=0; I<num_block_rows-1; I++){
          for(i=I*row_block_size, di=0; di<row_block_size; di++,i++){
              for(k=PB[i]; k<PE[i]; k++){
                  j = (*INDX)[k];
                  J = j / col_block_size;
                  if(block_count[J] == 0){
                      num_blocks ++;
                      block_count[J] ++;
                  }
              }
          }
          for(i=I*row_block_size, di=0; di<row_block_size; di++,i++){
              for(k=PB[i]; k<PE[i]; k++){
                  j = (*INDX)[k];
                  J = j / col_block_size;
                  block_count[J] = 0;
              }
          }   
      }
      for(i=I*row_block_size, di=0; di<tmp_num_rows_left; di++,i++){
          for(k=PB[i]; k<PE[i]; k++){
              j = (*INDX)[k];
              J = j / col_block_size;
              if(block_count[J] == 0){
                  num_blocks ++;
                  block_count[J] ++;
              }
          }
      }
      for(i=I*row_block_size, di=0; di<tmp_num_rows_left; di++,i++){
          for(k=PB[i]; k<PE[i]; k++){
              j = (*INDX)[k];
              J = j / col_block_size;
              block_count[J] = 0;
          }
      }   
    //End of Phase I

    *bnnz = num_blocks;
    *mb = num_block_rows;
    *kb = num_block_cols;
    
    bAi = (int*)aligned_malloc(sizeof(int)*num_blocks);
    bAj = (int*)aligned_malloc(sizeof(int)*num_blocks);
    bAx = (double*)aligned_malloc(sizeof(double)*num_blocks*block_size);
    memset (bAx, 0.0, sizeof(double)*num_blocks*block_size);
    double *blocks = (double*)aligned_malloc(row_block_size*col_block_size*num_block_cols*sizeof(double));
    memset (blocks, 0.0, sizeof(double)*row_block_size*col_block_size*num_block_cols);

    //Phase II: Copy all blocks.    
    int nnzIndex = 0;
    for(I=0; I<num_block_rows-1; I++){
        for(i=I*row_block_size, di=0; di<row_block_size; di++,i++){
            for(k=PB[i]; k<PE[i]; k++){
                j = (*INDX)[k];
                J = j / col_block_size;
                j0 = J * col_block_size;
                blocks[J*block_size + di*col_block_size + j-j0] = (*VAL)[k];
                block_count[J] ++;
            }
        }
        
        for(i=I*row_block_size, di=0; di<row_block_size; di++,i++){
            for(k=PB[i]; k<PE[i+1]; k++){
                j = (*INDX)[k];
                J = j / col_block_size;
                j0 = J * col_block_size;
                if(block_count[J] > 0){
                    memcpy(bAx+nnzIndex*block_size, blocks+J*block_size,sizeof(double)*block_size);
                    bAi[nnzIndex] = I;
                    bAj[nnzIndex] = J;
                    nnzIndex ++;
                    memset (blocks+J*block_size, 0.0, sizeof(double)*block_size); 
                    block_count[J] = 0;
                }
            }
        }
    }

    for(i=I*row_block_size, di=0; di<tmp_num_rows_left; di++,i++){
        for(k=PB[i]; k<PE[i]; k++){
            j = (*INDX)[k];
            J = j / col_block_size;
            j0 = J * col_block_size;
            blocks[J*block_size + di*col_block_size + j-j0] = (*VAL)[k];
            block_count[J] ++;
        }
    }

    for(i=I*row_block_size, di=0; di<tmp_num_rows_left; di++,i++){
        for(k=PB[i]; k<PE[i]; k++){
            j = (*INDX)[k];
            J = j / col_block_size;
            j0 = J * col_block_size;
            if(block_count[J] > 0){
                memcpy(bAx+nnzIndex*block_size, blocks+J*block_size,sizeof(double)*block_size);
                bAi[nnzIndex] = I;
                bAj[nnzIndex] = J;
                nnzIndex ++;
                memset (blocks+J*block_size, 0.0, sizeof(double)*block_size);
                block_count[J] = 0;
            }
        }  
    }
    //End of Phase II
    aligned_free(block_count);
    aligned_free(blocks);
    aligned_free(*INDX);
    aligned_free(*VAL);

    *INDX = bAi;             //copy the result to BCO_format matrix
    *JNDX = bAj;
    *VAL = bAx;
  }

