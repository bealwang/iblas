      module mod_lsbv_coo
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'COO'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************
      use mod_hash
      use mod_hash
      use representation_of_data
      use properties

      interface lsbv_coo
        module procedure ilsbv_coo
        module procedure slsbv_coo
        module procedure dlsbv_coo
        module procedure clsbv_coo
        module procedure zlsbv_coo
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_coo (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,base,ofs,nnz;
      char  diag,part;
      capsule* dummy;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       setup_hash(n,*ierr)
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i,*ierr)
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//do
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
       remove_hash(*ierr)
      }//function end  ilsbv_coo 
// **********************************************************************
// **********************************************************************
      void slsbv_coo (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,base,ofs,nnz;
      char  diag,part;
      capsule* dummy;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       setup_hash(n,*ierr)
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i,*ierr)
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0.0e0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//do
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0.0e0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
       remove_hash(*ierr)
      }//function end  slsbv_coo 
// **********************************************************************
// **********************************************************************
      void dlsbv_coo (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,base,ofs,nnz;
      char  diag,part;
      capsule* dummy;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       setup_hash(n,*ierr)
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i,*ierr)
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0.0d0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//do
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A[hash[i]->val_pos]!= 0.0d0 ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
       remove_hash(*ierr)
      }//function end  dlsbv_coo 
// **********************************************************************
// **********************************************************************
      void clsbv_coo (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,nnz;
      char  diag,part;
      capsule* dummy;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       setup_hash(n,*ierr)
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i,*ierr)
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A(hash[i]->val_pos)!= (0.0e0, 0.0e0) ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//do
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A(hash[i]->val_pos)!= (0.0e0, 0.0e0) ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
       remove_hash(*ierr)
      }//function end  clsbv_coo 
// **********************************************************************
// **********************************************************************
      void zlsbv_coo (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,nnz;
      char  diag,part;
      capsule* dummy;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='COO')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_infoa(mat->INFOA,'n',&nnz,*ierr);
      
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       setup_hash(n,*ierr)
      if (*ierr!=0) {//then
         return;
      }//end if
      for(i=0;i<n;i++){//donz
          new_capsule_main(mat->IA2[i]+ofs,mat->IA1[i]+ofs,i,*ierr)
         if (*ierr!=0) {//then
            return;
         }//end if
      }//end do
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A(hash[i]->val_pos)!= (0.0d0, 0.0d0) ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//do
            dummy=hash[i];
            while(dummy->pntr!=NULL){//while
               dummy=dummy->pntr;
               x[i]-=x[dummy->jndx] * mat->A[dummy->val_pos];
            }//end do
            if (diag!='U') {//then
               if(hash[i]->jndx==-1) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  if (mat->A(hash[i]->val_pos)!= (0.0d0, 0.0d0) ) {//then
                     x[i]=x[i]/mat->A[hash[i]->val_pos];
                  }else{/*else*/
                     *ierr = blas_error_singtria;
                     return;
                  }//end if
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
       remove_hash(*ierr)
      }//function end  zlsbv_coo 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_coo
