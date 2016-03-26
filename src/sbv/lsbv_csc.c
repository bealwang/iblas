      module mod_lsbv_csc
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'CSC'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************      
      use representation_of_data
      use properties

      interface lsbv_csc
        module procedure ilsbv_csc
        module procedure slsbv_csc
        module procedure dlsbv_csc
        module procedure clsbv_csc
        module procedure zlsbv_csc
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,base,ofs,pntr
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de=0;
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de=0;
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  ilsbv_csc 
// **********************************************************************
// **********************************************************************
      void slsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,base,ofs,pntr
      char  diag,part;
      real(KIND=sp)  :: de
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = 0.0e0 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0.0e0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = 0.0e0 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0.0e0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  slsbv_csc 
// **********************************************************************
// **********************************************************************
      void dlsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,base,ofs,pntr
      char  diag,part;
      real(KIND=dp)  :: de
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = 0.0d0 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0.0d0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = 0.0d0 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== 0.0d0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  dlsbv_csc 
// **********************************************************************
// **********************************************************************
      void clsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,pntr
      char  diag,part;
      complex(KIND=sp)  :: de
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = (0.0e0, 0.0e0) 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== (0.0e0, 0.0e0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = (0.0e0, 0.0e0) 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== (0.0e0, 0.0e0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  clsbv_csc 
// **********************************************************************
// **********************************************************************
      void zlsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,pntr
      char  diag,part;
      complex(KIND=dp)  :: de
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = (0.0d0, 0.0d0) 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== (0.0d0, 0.0d0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  pntr++;
               }//end do
            }else{/*else*/
               de = (0.0d0, 0.0d0) 
               pntr = mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  if(mat->IA1[pntr + ofs] + ofs==i) {//then
                     de=mat->A[pntr + ofs];
                  }else{/*else*/
                     x[i]-=mat->A[pntr + ofs] * x[mat->IA1[pntr + ofs ] + ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de== (0.0d0, 0.0d0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i] = x[i]/de
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zlsbv_csc 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_csc
