#include "INS_ROUTINER.h"
void dINS_entry(d_matrix * pmatrix, double val,int i,int j,int* istat)
{
  d_inelement* pelement;
  int nmb_element,ind;
  *istat=-1;
  if((i>pmatrix->DIM[1])||(j>pmatrix->DIM[2])) {//then
      *istat=blas_error_param;
      return;
    }
  else
    {
      nmb_element=new_d_element(pmatrix,istat);
      if(*istat!=0) return;
      ind=d_element_num(pmatrix,i,j,istat);
      if(*istat!=0) return;
      if(ind==0){// then  can't find the element mean this is new element
          pelement=daccess_element(nmb_element,pmatrix);
          if(*istat!=0) return;
          pelement->pntin.value=val;
          pelement->pntin.row_ind=i;
          pelement->pntin.col_ind=j;
        }else{
          pelement=daccess_element(ind,pmatrix);
          if(*istat!=0) return;
          pelement->pntin.value=pelement->pntin.value+val;
          dealloc_d_element(nmb_element,pmatrix,istat);
          if(*istat!=0) return;
        }
    }

}

void dINS_varbl_entr(d_matrix* vpmatrix,double val,int i,int j,int* istat)
{

  double** vall;
  int ii,jj,k,p,ind1,ind2,vall_ind1,vall_ind2;
  int m,n;
  int mm,nn;

  //determine the row of block entring
  ind1=-1;
  for(k=0;k<vpmatrix->DIM[3];k++)
    {
      ind1+=vpmatrix->sub_rows[k];
      if(ind1>=i) break;//the row length of block must bigger than the i
    }
  if(k<vpmatrix->DIM[3]) {//then
      ii=k;//block row index
      vall_ind1=(vpmatrix->sub_rows[k]-1)-(ind1-i) ;
    }else{
      *istat=blas_error_param;
      return;
    }//end if

  // determine the col of block entring
  ind2=-1;
  for(p=0;p<vpmatrix->DIM[4];p++)
    {
      ind2+=vpmatrix->sub_cols[p];
      if(ind2>=j) break;
    }//end do
  if(p<=vpmatrix->DIM[4]) {//then
      jj=p;//block col index
      vall_ind2=(vpmatrix->sub_cols[p]-1)-(ind2-j);
    }else{
      *istat=blas_error_param;
      return;
    }//end if

  mm=vpmatrix->sub_rows[ii];
  nn=vpmatrix->sub_cols[jj];
  vall=(double**)malloc2d(vpmatrix->sub_rows[ii],vpmatrix->sub_cols[jj],sizeof(double));
  
  #pragma omp parallel for num_threads(dtn(mm,MIN_ITERATOR_NUM)) private(n)
  for(m=0;m<mm;m++){
    for(n=0;n<nn;n++){
      vall[m][n]=0.0;
    }
  }

  vall[vall_ind1][vall_ind2]=val;
  dINS_varblock(vpmatrix,vall,vpmatrix->sub_rows[ii],vpmatrix->sub_cols[jj],ii,jj,istat);
  if(*istat!=0) return;

  aligned_free(vall);
  *istat=0;
}//end subroutine dINS_varbl_entr


void dINS_varblock(d_matrix* vpmatrix,double** val,int r_val,int c_val,int i,int j,int* istat)
{
  int nmb_element,ind;
  int m,n;
  //double** vv;
  d_inelement* pelement;
  int s_rows,s_cols,k;
  *istat=-1;
  s_rows=r_val;
  s_cols=c_val;
  //vv=(double**)malloc2d(s_rows,s_cols,sizeof(double));

  if((i>vpmatrix->DIM[3]||j>vpmatrix->DIM[4])||(s_rows!=vpmatrix->sub_rows[i])||(s_cols!=vpmatrix->sub_cols[j])){//then
      *istat=blas_error_param;
      return;
    }
  else
    {
      nmb_element=new_d_element(vpmatrix,istat);
      if(*istat!=0) return;
      ind=d_element_num(vpmatrix,i,j,istat);
      if(*istat!=0) return;
      if(ind==0)
        {//then
          pelement=daccess_element(nmb_element,vpmatrix);

          pelement->vblin.value=(double**)malloc2d (vpmatrix->sub_rows[i],vpmatrix->sub_cols[j],sizeof(double));

          #pragma omp parallel for num_threads(dtn(s_rows,MIN_ITERATOR_NUM)) private(n)
          for(m=0;m<s_rows;m++){
            for(n=0;n<s_cols;n++){
              pelement->vblin.value[m][n]=val[m][n];
            }
          }

          pelement->vblin.row_vblock_ind=i;
          pelement->vblin.col_vblock_ind=j;

          //update the row begin array and the row end array
          for(k=i;k<vpmatrix->DIM[3];k++)
            vpmatrix->trb[k]=vpmatrix->trb[k]+1;

          for(k=0;k<vpmatrix->DIM[3]-1;k++)
            vpmatrix->tre[k]=vpmatrix->trb[k+1];

          vpmatrix->tre[vpmatrix->DIM[3]-1]=nmb_element;
        }
      else
        {
          pelement=daccess_element(ind,vpmatrix);
          if(*istat!=0) return;

          for(m=0;m<s_rows;m++)
            for(n=0;n<s_cols;n++)
              pelement->vblin.value[m][n]+=val[m][n];

          dealloc_d_element(nmb_element,vpmatrix,istat);
          if(*istat!=0) return;
        }//end if
    }//end if
  *istat=0;
  return;
}//end subroutine dINS_varblock


void dINS_block (d_matrix* pmatrix,double** val,int r_val,int c_val,int i,int j,int* istat)
{
  int nmb_element,ind;
  //  double** vv;
  d_inelement* pelement;
  int s_rows,s_cols;
  int m,n;
  *istat=-1;
  s_rows=r_val;
  s_cols=c_val;
  //  vv=malloc2d (s_rows,s_cols,sizeof(**vv));

  if((i>pmatrix->DIM[3]||(j>pmatrix->DIM[4])||(s_rows!=pmatrix->DIM[5]||(s_cols!=pmatrix->DIM[6])))){//then
      *istat = blas_error_param;
      return;
    }else{
      nmb_element=new_d_element(pmatrix,istat);
      if(*istat!=0) return;
      ind=d_element_num(pmatrix,i,j,istat);
      if(*istat!=0) return;
      if(ind==0) {//then
          pelement=daccess_element(nmb_element,pmatrix);
          if(*istat!=0) return;
          pelement->blin.value=(double**)malloc2d (s_rows,s_cols,sizeof(double));

          for(m=0;m<s_rows;m++){
            for(n=0;n<s_cols;n++){
              pelement->blin.value[m][n]=val[m][n];
            }
          }

          pelement->blin.row_block_ind=i;
          pelement->blin.col_block_ind=j;
        }else{
          pelement=daccess_element(ind,pmatrix);
          if(*istat!=0) return;

          #pragma omp parallel for num_threads(dtn(s_rows,MIN_ITERATOR_NUM)) private(n)
          for(m=0;m<s_rows;m++){
            for(n=0;n<s_cols;n++)
              {
                //                vv[m][n]=vv[m][n]+val[m][n];//??? what does it use for?
                pelement->blin.value[m][n]+=val[m][n];
              }
            }
          dealloc_d_element (nmb_element,pmatrix,istat);
          if(*istat!=0) return;
        }//end if
    }//end if

  //  free(vv);
  *istat=0;
}//end subroutine dINS_block


void dINS_bl_entr (d_matrix* pmatrix,double val,int i,int j,int* istat)
{
  double** vall;
  int ii,jj;
  int m,n;
  *istat=-1;
  ii=i/(pmatrix->DIM[5]);
  jj=j/(pmatrix->DIM[6]);


  vall=(double**)malloc2d(pmatrix->DIM[5],pmatrix->DIM[6],sizeof(double));

  for(m=0;m<pmatrix->DIM[5];m++){
    for(n=0;m<pmatrix->DIM[6];n++){
      vall[m][n]=0;
    }
  }

  vall[i-ii*pmatrix->DIM[5]][j-jj*pmatrix->DIM[6]]=val;
  dINS_block (pmatrix,vall,pmatrix->DIM[5],pmatrix->DIM[6],ii,jj,istat);

  if(*istat!=0) return;
  aligned_free(vall);
  *istat=0;
}//end subroutine dINS_bl_entr


dsp_linknode* duscr_blockend(d_matrix* A,int prpty,int* istat)
{
  int m,n,bnnz,ind,kb,lb,dummy,k,mb,i,j;
  int *bindx,*bjndx;
  double* val;
  int nmb_block;
  d_matrix* pmatrix;
  d_inelement* pelement;
  dsp_linknode* dsp_l;
  *istat=-1;
  pmatrix=daccess_matrix(A);
  if(pmatrix==NULL) return NULL;
  lb=pmatrix->DIM[5];//the row size of the block
  bnnz=pmatrix->d_element_start->number;//the number of total blocks
  mb=pmatrix->DIM[3];//the number of block in one row
  kb=pmatrix->DIM[4];//the number of block in one colume
  m=mb*lb;//the width of the matrix
  n=kb*lb;//the height of the matrix
  ind=0;
  dummy=bnnz*lb*lb;//the size of space which can stored bnnz blocks with echo block in lb*lb size
  val=(double*)aligned_malloc(dummy*sizeof(double));
  bindx=(int*)aligned_malloc(bnnz*sizeof(int));
  bjndx=(int*)aligned_malloc(bnnz*sizeof(int));

  nmb_block=bnnz;

  for(i=0;i<bnnz;i++)
    {
      //  k=1;
      pelement=daccess_element(nmb_block-i,pmatrix);//num from 1 to bnnz

      for(j=0;j<lb;j++){
          for(k=0;k<lb;k++)
            {
              ind++;
              val[ind]=pelement->blin.value[k][j];//颠倒了？这里是列优先
            }//end do
        }//end do
      bindx[i]=pelement->blin.row_block_ind;
      bjndx[i]=pelement->blin.col_block_ind;
    }//end do

  // RELEASING
  d_dealloc (A,istat);
  if (*istat!=0) return NULL;
  // CREATE A MATRIX IN BCO FORMAT
  *istat=-1 ;//needed to create copy of data
  dsp_l=duscr_bco(m,n,val,dummy,bindx,bnnz,bjndx,bnnz,bnnz,mb,kb,lb,prpty,istat);

  if (*istat!=0) return NULL;

  aligned_free(val);
  aligned_free(bindx);
  aligned_free(bjndx);
  return dsp_l;
}//end subroutine  duscr_blockend

dsp_linknode* duscr_varend (d_matrix* A,int prpty,int *istat)
{
  int  m,n,ind,kb,mb;
  int *bindx,*indx,*rpntr,*cpntr,*bpntrb,*bpntre;
  double* val;
  int size_val,val_ind,indx_ind,bindx_ind,ii,jj,i,j;
  d_matrix* pmatrix;
  d_inelement* pelement;
  dsp_linknode* dsp_l;
  *istat=-1;
  pmatrix=daccess_matrix(A);

  mb=pmatrix->DIM[3];
  kb=pmatrix->DIM[4];
  // determine  size of val,m,n
  size_val=0;
  m=0;
  n=0;


  for(i=0;i<pmatrix->DIM[3];i++)
    {
      for(j=0;j<pmatrix->DIM[4];j++)
        {
          ind=d_element_num (pmatrix,i,j,istat);
          if(ind!=0) {//then
              pelement=daccess_element(ind,pmatrix);
              if(pelement!=NULL)
                size_val+=pmatrix->sub_rows[i]*pmatrix->sub_cols[j];
            }//end if
        }//end do
    }//end do


  for(i=0;i<pmatrix->DIM[3];i++){
      m=m+pmatrix->sub_rows[i];
    }//total rows
  for(j=0;j<pmatrix->DIM[4];j++){
      n=n+pmatrix->sub_cols[j];
    }//total clos
  val=(double*)aligned_malloc(sizeof(double)*size_val);
  indx=(int*)aligned_malloc(sizeof(int)*(pmatrix->d_element_start->number+1));
  bindx=(int*)aligned_malloc(sizeof(int)*(pmatrix->d_element_start->number));

  for(i=0;i<size_val;i++)
    val[i]=0.0;

  //fill val ,indx and bindx
  val_ind=0;
  indx_ind=1;
  bindx_ind=0;
  indx[0]=0;
  for(i=0;i<pmatrix->DIM[3];i++){//for
      for(j=0;j<pmatrix->DIM[4];j++){//for
          ind=d_element_num(pmatrix,i,j,istat);
          if (*istat!=0) return NULL;
          if(ind!=0) {//then
              pelement=daccess_element(ind,pmatrix);
              if (pelement==NULL) return NULL;
              for(jj=0;jj<pmatrix->sub_cols[j];jj++){//for
                  for(ii=0;ii<pmatrix->sub_rows[i];ii++){//for
                      val[val_ind]=pelement->vblin.value[ii][jj];
                      val_ind++;
                    }//end do
                }//end do
              bindx[bindx_ind]=j;
              bindx_ind++;
              indx[indx_ind]=indx[indx_ind-1]+pmatrix->sub_rows[i]*pmatrix->sub_cols[j];
              indx_ind++;//块的index范围
            }//end if
        }//end do
    }//end do




  //fill rpntr, cpntr,bpntrb,bpntre
  rpntr=(int*)aligned_malloc(sizeof(int)*(pmatrix->DIM[3]+1));
  cpntr=(int*)aligned_malloc(sizeof(int)*(pmatrix->DIM[4]+1));
  bpntrb=(int*)aligned_malloc(sizeof(int)*(pmatrix->DIM[3]));
  bpntre=(int*)aligned_malloc(sizeof(int)*(pmatrix->DIM[3]));

  rpntr[0]=0;
  for(i=1;i<pmatrix->DIM[3]+1;i++){
      rpntr[i]=rpntr[i-1]+pmatrix->sub_rows[i-1];
    }//end do

  cpntr[0]=0;
  for(j=1;j<pmatrix->DIM[4]+1;j++){
      cpntr[j]=cpntr[j-1]+pmatrix->sub_cols[j-1];
    }//end do

  for(i=0;i<pmatrix->DIM[3];i++){
      bpntrb[i]=pmatrix->trb[i];
      bpntre[i]=pmatrix->tre[i];
    }//end do



  //CREATING MATRIX IN VBR FORMAT
  *istat=-1; //needed to create copy of data
  dsp_l=duscr_vbr(m,n,val,size_val,indx,(pmatrix->d_element_start->number+1),bindx,(pmatrix->d_element_start->number),rpntr,(pmatrix->DIM[3]+1),cpntr,(pmatrix->DIM[4]+1),bpntrb,(pmatrix->DIM[3]),bpntre,(pmatrix->DIM[3]),mb,kb,prpty,istat);
  if (*istat!=0) return NULL;

  aligned_free(val);
  aligned_free(bindx);
  aligned_free(indx);
  aligned_free(rpntr);
  aligned_free(cpntr);
  aligned_free(bpntrb);
  aligned_free(bpntre);

  //RELEASING
  d_dealloc (pmatrix,istat);
 if (*istat!=0) return NULL;

  *istat=0;
  return dsp_l;
}//end subroutine  duscr_varend

dsp_linknode* duscr_normend(d_matrix* A,int prpty,int *istat)
{
  int  m,n,nnz;
  int *indx,*jndx;
  double* val;
  int nmb_element,i;
  d_matrix *pmatrix;
  d_inelement*pelement;
  dsp_linknode* dsp_l;
  *istat=-1;
  pmatrix=daccess_matrix(A);

  m=pmatrix->DIM[1];//nb_of_rows
  n=pmatrix->DIM[2];//nb_of_cols
  nnz=pmatrix->d_element_start->number;


  val=(double*)aligned_malloc(sizeof(double)*nnz);
  indx=(int*)aligned_malloc(sizeof(int)*nnz);
  jndx=(int*)aligned_malloc(sizeof(int)*nnz);
  nmb_element=nnz;

  for(i=0;i<nnz;i++)
    {
      pelement=daccess_element(nmb_element-i,pmatrix);
      if (pelement==NULL) return NULL;
      val[i]=pelement->pntin.value;
      indx[i]=pelement->pntin.row_ind;
      jndx[i]=pelement->pntin.col_ind;
    }//end do

  d_dealloc (A,istat);
  if (*istat!=0) return NULL;

  //CREATING A MATRIX IN COO FORMAT
  *istat=-1; //needed to create copy of data
  dsp_l=duscr_coo(m,n,val,nnz,indx,jndx,nnz,prpty,istat);
  if (*istat!=0) return NULL;
  aligned_free(val);
  aligned_free(indx);
  aligned_free(jndx);


  *istat=0;
  return dsp_l;
}//end  subroutine duscr_normend
