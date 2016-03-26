      module mod_lsbv_bdi
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BDI'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************
            
      interface lsbv_bdi
        module procedure ilsbv_bdi
        module procedure slsbv_bdi
        module procedure dlsbv_bdi
        module procedure clsbv_bdi
        module procedure zlsbv_bdi
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_bdi (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,j,n,mm,nn,blda,nbdiag,dd,nn_sq;
      char diag,part,store;
      int* y;
      *ierr=-1;
      n=n_x;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&blda,*ierr);
      
       get_infoa(mat->INFOA,'g',&nbdiag,*ierr);
      
      *ierr=-1;
      if ((mat->FIDA!='BDI')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0 
      if (part=='L') {//then
         for(i=blda-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_left_lower(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*(dd-1)+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=0;i<blda;i++){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_right_upper(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*dd+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
      }//function end  ilsbv_bdi 
// **********************************************************************
// **********************************************************************
      void slsbv_bdi (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,j,n,mm,nn,blda,nbdiag,dd,nn_sq;
      char diag,part,store;
      float* y;
      *ierr=-1;
      n=n_x;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&blda,*ierr);
      
       get_infoa(mat->INFOA,'g',&nbdiag,*ierr);
      
      *ierr=-1;
      if ((mat->FIDA!='BDI')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0e0 
      if (part=='L') {//then
         for(i=blda-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_left_lower(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*(dd-1)+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=0;i<blda;i++){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_right_upper(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*dd+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
      }//function end  slsbv_bdi 
// **********************************************************************
// **********************************************************************
      void dlsbv_bdi (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,j,n,mm,nn,blda,nbdiag,dd,nn_sq;
      char diag,part,store;
      double* y;
      *ierr=-1;
      n=n_x;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&blda,*ierr);
      
       get_infoa(mat->INFOA,'g',&nbdiag,*ierr);
      
      *ierr=-1;
      if ((mat->FIDA!='BDI')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0d0 
      if (part=='L') {//then
         for(i=blda-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_left_lower(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*(dd-1)+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=0;i<blda;i++){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_right_upper(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*dd+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
      }//function end  dlsbv_bdi 
// **********************************************************************
// **********************************************************************
      void clsbv_bdi (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,mm,nn,blda,nbdiag,dd,nn_sq;
      char diag,part,store;
      complex_f* y;
      *ierr=-1;
      n=n_x;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&blda,*ierr);
      
       get_infoa(mat->INFOA,'g',&nbdiag,*ierr);
      
      *ierr=-1;
      if ((mat->FIDA!='BDI')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0e0,0.0e0) 
      if (part=='L') {//then
         for(i=blda-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_left_lower(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*(dd-1)+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=0;i<blda;i++){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_right_upper(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*dd+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
      }//function end  clsbv_bdi 
// **********************************************************************
// **********************************************************************
      void zlsbv_bdi (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,mm,nn,blda,nbdiag,dd,nn_sq;
      char diag,part,store;
      complex_d* y;
      *ierr=-1;
      n=n_x;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&blda,*ierr);
      
       get_infoa(mat->INFOA,'g',&nbdiag,*ierr);
      
      *ierr=-1;
      if ((mat->FIDA!='BDI')||(mat->M!=n)||(mat->K!=n)||(mm!=nn)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0d0,0.0d0) 
      if (part=='L') {//then
         for(i=blda-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((-mat->IA1[j]<blda-i+1)&&(mat->IA1[j]<0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_left_lower(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*(dd-1)+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=0;i<blda;i++){//for
            if (diag=='U') {//then
               for(j=0;j<nb;j++){//end fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                  x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
            }else{/*else*/
               dd=-1;
               for(j=0;j<nb;j++){//end fordiag
                  if (mat->IA1[j]==0) {//then
                     dd=j;
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
                      block_T_mult_vec(mat->A[(blda*j-mat->IA1[j]+i-1)*nn_sq+1][(blda*(j-1)-mat->IA1[j]+i)*nn_sq],x[(i-mat->IA1[j]-1)*nn+1][(i-mat->IA1[j])*nn],nn,y,nn,store,*ierr);
                     x[i*nn+1][i*nn]-=y;
                  }//end if
               }//end do
               if (dd!=-1) {//then 
                   invert_T_right_upper(mat->A[(blda*dd+i-1)*nn_sq+1][(blda*dd+i)*nn_sq],x[i*nn+1][i*nn],nn,store,*ierr);
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr = blas_error_memdeloc;
            return;
         }//end if
      }//end if
      }//function end  zlsbv_bdi 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_bdi
