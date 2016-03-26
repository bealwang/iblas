#include <assert.h>
#include "comm_tools.h"
#include "INSERTING.h"
#include "types.h"
#include "properties.h"
#include "link.h"
#include "conv/conv_tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*计算循环需要的线程数量
根据循环迭代次数和CPU核数及一个线程最少需要的循环迭代次数
来计算出需要的线程数量，计算出的最大线程数量不超过CPU核数
@param int n  -循环迭代次数
@param int min_n  -单个线程需要的最少迭代次数
@return int   -线程数量
*/
int dtn(int n,int min_n){
  const int g_ncore = omp_get_num_procs();
  int max_tn = n/min_n;
  int tn = max_tn>3*g_ncore?3*g_ncore:max_tn;
  if(tn < 1){
    tn = 1;
  }
  return tn;
}

int count (double* a, int n_a, int e)
{
  int i=0;
  int counter=0;
      for(i=0;i<n_a;i++)
        if(a[i]==e) counter++;
      return counter;
}

// void collect_features(dsp_linknode* spmat,SpFeature_CPU *spf,int *ierr){
//   DSPMAT*  dsp_data;
//   int num_rows,num_cols,num_nonzeros,i,j,max_RD,temp;
//   double squre_sum = 0.0,aver_RD,var_RD;
//   dsp_data = accessdata_dsp(spmat,ierr);
//   if (*ierr!=0) {
//     *ierr = blas_error_param;
//     printf("2\n");
//     return;
//   }
//   num_rows = dsp_data->spc.M;
//   num_cols = dsp_data->spc.N;
//   num_nonzeros = dsp_data->spc.NNZ;
  
//   (*spf).M = num_rows;
//   (*spf).N = num_cols;
//   (*spf).NNZ = num_nonzeros;

//   aver_RD = (double) num_nonzeros/num_rows;
//   (*spf).avg_RD = aver_RD;

//   usconv_coo2csr(spmat,ierr);
//   if(0 != *ierr){
//     return;
//   }
//   dsp_data = accessdata_dsp(spmat,ierr);
//   int* nzs_per_row = (int*)aligned_malloc(num_rows*sizeof(int));
//   memset(nzs_per_row,0,num_rows*sizeof(int));
//   temp = max_RD = 0;
//   for(i=0; i<dsp_data->M; i++){
//     temp = dsp_data->PE[i] - dsp_data->PB[i];
//     if(max_RD < temp){
//       max_RD = temp;
//     }
//     nzs_per_row[i] = temp;
//     squre_sum += (temp-aver_RD) * (temp-aver_RD);
//     temp = 0;
//   }

//   (*spf).max_RD = max_RD;
//   (*spf).var_RD = squre_sum/num_rows;

//   usconv_csr2coo(spmat,ierr);
//   if(0 != *ierr){
//     return;
//   }
//   dsp_data = accessdata_dsp(spmat,ierr);
//   /*-------------------------计算幂律特征---------------------------*/
//     //Calculate the pow-law parameters
//     double *nzs_distribution = (double*)aligned_malloc((num_cols+1)*sizeof(double));
//     memset(nzs_distribution,0.0,(num_cols+1)*sizeof(double));
//     int number;
//     for(i = 0; i < num_rows; i++)
//     {
//         number = nzs_per_row[i];
//         nzs_distribution[number] += 1;
//     }
//     aligned_free(nzs_per_row);
    
//     int total_count = 0;
//     int peak_pos = 0;
//     double total_sum = 0, peak_ratio = 0;
//     for(i = 1; i <= num_cols; i++)
//     {
//         if ( nzs_distribution[i] != 0.0)
//         {
//             total_count ++;
//             total_sum += nzs_distribution[i];
//             if ( peak_ratio < nzs_distribution[i] )
//             {
//               peak_ratio = nzs_distribution[i];
//               peak_pos = total_count;
//             }
//         }
//     }
//     double *final_Pks = (double*)aligned_malloc((total_count+1)*sizeof(double));
//     int *final_ks = (int*)aligned_malloc((total_count+1)*sizeof(int));
//     memset(final_Pks,0.0,(total_count+1)*sizeof(double));
//     memset(final_ks,0,(total_count+1)*sizeof(int));
//     for (i = 1, j = 1; i <= num_cols; i++)
//     {
//         if ( nzs_distribution[i] != 0.0)
//         {
//             final_Pks[j] = nzs_distribution[i];
//             final_ks[j] = i;
//             j ++;
//         }
//     }
//     aligned_free(nzs_distribution);

//     int count = 0;
//     double aver_x = 0.0, aver_y = 0.0;
//     for(i = peak_pos; i <= total_count; i++)
//     {
//         count ++;
//         aver_x += log10(final_ks[i]);
//         aver_y += log10(final_Pks[i]);

//         aver_x += 1;
//         aver_y += 1;
//     }
//     assert ( count == (total_count - peak_pos + 1) );
//     aver_x /= count;
//     aver_y /= count;
//     //printf("aver_x:%lf\naver_y:%lf\n",aver_x,aver_y);
//     double a_up = 0.0, a_down = 0.0;
//     double a = 0.0;
//     for(i = peak_pos; i <= total_count; i++)
//     {
//         a_up += (log10(final_ks[i]) - aver_x) * (log10(final_Pks[i]) - aver_y);
//         a_down += (log10(final_ks[i]) - aver_x) * (log10(final_ks[i]) - aver_x);  //log10(final_ks[3] == aver_x)!!!!
//         // a_up += 1;
//         // a_down += 1;
//     }
//     aligned_free(final_ks);
//     aligned_free(final_Pks);
//     //printf("a_down:%lf\n",a_down);
//     a = a_up / a_down;
//     (*spf).R = 0.0-a;
//   /*----------------------------------计算幂律特征------------------------*/

//   /*-----------------DIA--------------------*/
//   size_t complete_ndiags = 0;
//   int unmarked = -1; // works for both signed and unsigned
//   int ii,jj,map_index;
//   int* diag_map = (int*)aligned_malloc(sizeof(int)*(num_rows+num_cols));     //mark the diags
//   int* diag_map_2 = (int*)aligned_malloc(sizeof(int)*(num_cols+num_rows));
//   memset(diag_map,unmarked,sizeof(int)*(num_rows+num_cols));
//   memset(diag_map_2,0,sizeof(int)*(num_rows+num_cols));
//   for(i=0;i<num_nonzeros;i++){
//     ii = dsp_data->IA1[i];
//     jj = dsp_data->IA2[i];
//     map_index = num_rows-ii+jj;            //used to find the same diag
//     if(diag_map[map_index] == unmarked){        
//       diag_map[map_index] = complete_ndiags;
//       complete_ndiags++;                 //number of diags
//     }
//     diag_map_2[map_index] ++;
//   }
//   aligned_free(diag_map);

//   int j_ndiags = 0;
//   double temp_ratio;
//   int NTdiags = 0;
//   double *array_ndiags = (double*)aligned_malloc(10*sizeof(double));
//   memset(array_ndiags,0,sizeof(double)*10);
//   for (i=0; i<num_cols+num_rows; i++)
//   {
//       if (diag_map_2[i] != 0)
//       {
//           j_ndiags ++;
//           temp_ratio = (double)diag_map_2[i] / num_rows;
//           if (temp_ratio < 0.1 )
//             array_ndiags[0] ++;
//           else if (temp_ratio < 0.2 )
//             array_ndiags[1] ++;
//           else if (temp_ratio < 0.3 )
//             array_ndiags[2] ++;
//           else if (temp_ratio < 0.4 )
//             array_ndiags[3] ++;
//           else if (temp_ratio < 0.5 )
//             array_ndiags[4] ++;
//           else if (temp_ratio < 0.6 )
//             array_ndiags[5] ++;
//           else if (temp_ratio < 0.7 )
//             array_ndiags[6] ++;
//           else if (temp_ratio < 0.8 )
//             array_ndiags[7] ++;
//           else if (temp_ratio < 0.9 )
//             array_ndiags[8] ++;
//           else if (temp_ratio <= 1.0 )
//             array_ndiags[9] ++;

//           if (temp_ratio >= 0.6 )
//               NTdiags ++;
//       }
//   }
//   assert( j_ndiags == complete_ndiags);
//   aligned_free(diag_map_2);

//   for ( i=0; i<10; i++)
//   {
//     array_ndiags[i] /= complete_ndiags;
//   }
//   double NTdiags_ratio = (double)NTdiags/complete_ndiags;
//   double ER_DIA = (double)num_nonzeros / (complete_ndiags * num_rows);

//   aligned_free(array_ndiags);

//   (*spf).Ndiags = complete_ndiags;
//   (*spf).NTdiags_ratio = NTdiags_ratio;
//   (*spf).er_dia = ER_DIA;
//   /*--------------------------DIA------------------*/
// }

void showFeature(SpFeature_CPU spf){
  printf("---------  features  ---------------\n");
  printf("行数：  %d\n",spf.M);
  printf("列数：  %d\n",spf.N);
  printf("非零元个数：  %d\n",spf.NNZ);
  printf("行非零元个数最大值：  %d\n",spf.max_RD);
  printf("行非零元个数平均值：  %lf\n",spf.avg_RD);
  printf("行非零元波动方差：  %lf\n",spf.var_RD);
  printf("对角线条数：%d\n",spf.Ndiags);
  printf("真对角线比例:  %lf\n",spf.NTdiags_ratio);
  printf("DIA格式非零元比例:  %lf\n",spf.er_dia);
  printf("幂律分布因子：  %lf\n",spf.R);
  printf("-------  end feature  --------------\n");
}

void* aligned_malloc(size_t size){
  void* non_aligned;
  void** aligned;
  non_aligned = (void*)malloc((size_t)size+(MEM_ALIGNMENT-1)+sizeof(void*));
  assert(non_aligned != NULL);
  aligned = (void**)( (size_t)(((size_t)(non_aligned)+MEM_ALIGNMENT+sizeof(void*)) & ~(MEM_ALIGNMENT-1)));
  *(aligned-1) = non_aligned;

  return (void*)aligned;

  // void* a;
  // a = (void*)malloc(size);
  // assert(a != NULL);
  // return a;
}

void aligned_free(void* p){
  assert(p != NULL);
  free(*((void **)p-1));
  // free(p);
}
void** malloc2d(int rows, int cols, int size)
{
  int j;
  int rowSize = cols * size;
  int indexSize = rows * sizeof(void *);
  void **a = (void **)aligned_malloc(indexSize + rows* rowSize);
  assert(a != NULL);
  char *dataStart = (char *) a + indexSize;
  for(j = 0; j < rows; j++){
  a[j] = dataStart + j * rowSize;
  }
  
  return a;
}

void initarray(int **array,int rows,int cols){
  int i,j;
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      array[i][j] = 0;
    }
  }
}
int minval(int *a, int n_a)
{
  int x,i;
  if(n_a<=0) return -1;//
  x=a[0];
  for(i=1;i<n_a;i++){
    if(a[i]<x)
      x=a[i];
  }
  return x;
}
int maxval (int *a, int n_a)
{
  int x,i;
  if(n_a<=0) return -1;
  x=a[0];
  for(i=1;i<n_a;i++)
    if(a[i]>x)
      x=a[i];
  return x;
}

void dump_matrix (void* matrix)
{
  d_matrix* p=matrix;
  d_element* e=p->d_element_start;
  int i,j;
  double** o;
  switch(p->format)
    {
    case 'n':
      printf("dump begin===================================================\n");
      printf("format:normal\n");
      printf("size:%d X %d\n\n",p->DIM[1],p->DIM[2]);
      o=(double**)malloc2d(p->DIM[1],p->DIM[2],sizeof(double));
      #pragma omp parallel for num_threads(dtn(p->DIM[1],MIN_ITERATOR_NUM)) private(j)
      for(i=0;i<p->DIM[1];i++){
        for(j=0;j<p->DIM[2];j++){
          o[i][j]=0;
        }
      }
      
      for(i=0;e!=NULL;i++)
        {
          o[e->contents.pntin.row_ind][e->contents.pntin.col_ind]=e->contents.pntin.value;
          e=e->pntr;
        }
      for(i=0;i<p->DIM[1];i++)
        {
          for(j=0;j<p->DIM[2];j++){
            printf("%2.1f ",o[i][j]);
          }
          printf("\n");
        }
      printf("dump end===================================================\n\n");
      aligned_free(o);
      break;
    case 'b':
      printf("dump begin===================================================\n\n");
      printf("format:block\n");
      printf("size:block(%d,%d)in %d X %d\n",p->DIM[5],p->DIM[6],p->DIM[3],p->DIM[4]);

      o=(double**)malloc2d(p->DIM[3]*p->DIM[5],p->DIM[4]*p->DIM[6],sizeof(double));
      #pragma omp parallel for num_threads(dtn(p->DIM[3]*p->DIM[5],MIN_ITERATOR_NUM)) private(j)
      for(i=0;i<p->DIM[3]*p->DIM[5];i++)
        for(j=0;j<p->DIM[4]*p->DIM[6];j++)
          o[i][j]=0;

      while(e)
        {
          for(i=0;i<p->DIM[5];i++)
            for(j=0;j<p->DIM[6];j++)
              o[e->contents.blin.row_block_ind*p->DIM[5]+i][e->contents.blin.col_block_ind*p->DIM[6]+j]=e->contents.blin.value[i][j];
        }
      for(i=0;i<p->DIM[3]*p->DIM[5];i++)
       {
          for(j=0;j<p->DIM[4]*p->DIM[6];j++)
            printf("%2.1f ",o[i][j]);
          printf("\n");
        }
      printf("dump begin===================================================\n\n");
      aligned_free(o);
      break;
    case 'v':
      printf("vblock\n");
      break;
    default:
      break;
    }
}
