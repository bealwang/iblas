      module mod_lsbv_bsr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'BSR'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************      
                  
      interface lsbv_bsr
        module procedure ilsbv_bsr
        module procedure slsbv_bsr
        module procedure dlsbv_bsr
        module procedure clsbv_bsr
        module procedure zlsbv_bsr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
int j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq;
      char diag,part,store;
      int* y;
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
      bofs=-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0 
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs]+ ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function ilsbv_bsr 
// **********************************************************************
// **********************************************************************
      void slsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq;
      char diag,part,store;
      float* y;
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
      bofs=-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0e0 
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs]+ ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function slsbv_bsr 
// **********************************************************************
// **********************************************************************
      void dlsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
      int j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq;
      char diag,part,store;
      double* y;
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
      bofs=-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0d0 
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs]+ ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function dlsbv_bsr 
// **********************************************************************
// **********************************************************************
      void clsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      int j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq;
      char diag,part,store;
      complex_f* y;
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
      bofs=-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0e0, 0.0e0) 
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs]+ ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function clsbv_bsr 
// **********************************************************************
// **********************************************************************
      void zlsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      int j,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq;
      char diag,part,store;
      complex_d* y;
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
      bofs=-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0d0, 0.0d0) 
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs]+ ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                   block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
while((pntr<mat->pe[j])&&(mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  dd=pntr;
               }//end if
               if(dd==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[(dd + bofs)*nn_sq+1][(dd + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr + ofs] + ofs!=j) {//then
                      block_T_mult_vec(mat->A[(pntr + bofs)*nn_sq+1][(pntr + bofs +1)*nn_sq],x[j*nn+1][j*nn],nn,y,nn,store,*ierr);  
      x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn]-=y;
                  }//end if
                  pntr++;
               }//end do
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function zlsbv_bsr 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_bsr
