      module mod_lsbv_bco
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'BCO'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************
                        interface lsbv_bco
        module procedure ilsbv_bco
        module procedure slsbv_bco
        module procedure dlsbv_bco
        module procedure clsbv_bco
        module procedure zlsbv_bco
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_bco (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,ofs,base,nb,mb,nnz;
      int  mm,nn,nn_sq;
      char diag,part,store;
      capsule* dummy;
      int* y; //extra stor.//
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      if ((mat->FIDA!='BCO')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0 
       setup_hash(nb,*ierr);
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i*nn_sq+1,*ierr);
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=nb-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
         for(i=0;i<n;i++){//forb
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
       remove_hash(*ierr)
      }//function end  ilsbv_bco 
// **********************************************************************
// **********************************************************************
      void slsbv_bco (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,ofs,base,nb,mb,nnz;
      int  mm,nn,nn_sq;
      char diag,part,store;
      capsule* dummy;
      float* y; //extra stor.//
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      if ((mat->FIDA!='BCO')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0e0 
       setup_hash(nb,*ierr);
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i*nn_sq+1,*ierr);
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=nb-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
         for(i=0;i<n;i++){//forb
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
       remove_hash(*ierr)
      }//function end  slsbv_bco 
// **********************************************************************
// **********************************************************************
      void dlsbv_bco (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,ofs,base,nb,mb,nnz;
      int  mm,nn,nn_sq;
      char diag,part,store;
      capsule* dummy;
      double* y; //extra stor.//
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      if ((mat->FIDA!='BCO')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = 0.0d0 
       setup_hash(nb,*ierr);
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i*nn_sq+1,*ierr);
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=nb-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
         for(i=0;i<n;i++){//forb
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
       remove_hash(*ierr)
      }//function end  dlsbv_bco 
// **********************************************************************
// **********************************************************************
      void clsbv_bco (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,ofs,base,nb,mb,nnz;
      int  mm,nn,nn_sq;
      char diag,part,store;
      capsule* dummy;
      complex_f* y; //extra stor.//
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      if ((mat->FIDA!='BCO')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0e0, 0.0e0) 
       setup_hash(nb,*ierr);
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i*nn_sq+1,*ierr);
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=nb-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
         for(i=0;i<n;i++){//forb
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
       remove_hash(*ierr)
      }//function end  clsbv_bco 
// **********************************************************************
// **********************************************************************
      void zlsbv_bco (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,ofs,base,nb,mb,nnz;
      int  mm,nn,nn_sq;
      char diag,part,store;
      capsule* dummy;
      complex_d* y; //extra stor.//
      *ierr=-1;
      n=n_x;
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'d',&mm,*ierr);
      
       get_infoa(mat->INFOA,'e',&nn,*ierr);
      
       get_infoa(mat->INFOA,'f',&mb,*ierr);
      
       get_infoa(mat->INFOA,'g',&nb,*ierr);
      
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'f',&store,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      if ((mat->FIDA!='BCO')||(mat->M!=n)||(mat->K!=n)|| (mm!=nn)||(nb!=mb)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      nn_sq=nn*nn;
      y=malloc(nn*sizeof(y*));
      
      y = (0.0d0, 0.0d0) 
       setup_hash(nb,*ierr);
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i*nn_sq+1,*ierr);
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=nb-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_left_lower(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
         for(i=0;i<n;i++){//forb
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
                block_T_mult_vec(mat->A[dummy->val_pos][dummy->val_pos+nn_sq-1],x[(dummy->jndx-1)*nn+1][(dummy->jndx)*nn],nn,y,nn,store,*ierr);
               x[i*nn+1][i*nn]-=y;
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                   invert_T_right_upper(mat->A[hash[i]->val_pos][hash[i]->val_pos+nn_sq-1],x[i*nn+1][i*nn],nn,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
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
       remove_hash(*ierr)
      }//function end  zlsbv_bco 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_bco
