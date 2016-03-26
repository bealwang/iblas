      module mod_rsbv_bsr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BSR'-STORAGE
//                   rsbv=Right Solve By Vector
// **********************************************************************
                  
      interface rsbv_bsr
        module procedure irsbv_bsr
        module procedure srsbv_bsr
        module procedure drsbv_bsr
        module procedure crsbv_bsr
        module procedure zrsbv_bsr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void irsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
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
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)||(nb!=mb)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y=0 
      if (part=='L') {//then
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_left_lower(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_right_upper(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function irsbv_bsr 
// **********************************************************************
// **********************************************************************
      void srsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
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
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)||(nb!=mb)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y=0.0e0 
      if (part=='L') {//then
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_left_lower(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_right_upper(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function srsbv_bsr 
// **********************************************************************
// **********************************************************************
      void drsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
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
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)||(nb!=mb)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y=0.0d0 
      if (part=='L') {//then
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_left_lower(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_right_upper(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function drsbv_bsr 
// **********************************************************************
// **********************************************************************
      void crsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
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
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)||(nb!=mb)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y=(0.0e0,0.0e0) 
      if (part=='L') {//then
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_left_lower(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_right_upper(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function crsbv_bsr 
// **********************************************************************
// **********************************************************************
      void zrsbv_bsr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,pntr,ofs,bofs,mm,nn,mb,nb,dd,nn_sq
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
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if ((mat->FIDA!='BSR')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)||(nb!=mb)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y=(0.0d0,0.0d0) 
      if (part=='L') {//then
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_left_lower(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                   block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                  x[i*nn+1][i*nn]-=y;
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                      block_mult_vec(mat->A[(pntr+bofs)*nn_sq+1][(pntr+bofs+1)*nn_sq],x[(mat->IA1[pntr+ofs]+bofs)*nn+1][(mat->IA1[pntr+ofs]+bofs+1)*nn],nn,y,nn,store,*ierr); 
                     x[i*nn+1][i*nn]-=y;
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_right_upper(mat->A[(dd+bofs)*nn_sq+1][(dd+bofs+1)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function zrsbv_bsr 
// **********************************************************************
// **********************************************************************
      end module mod_rsbv_bsr
