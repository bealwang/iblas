      module mod_rsbv_dia
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'DIA'-STORAGE
//                   rsbv=Right Solve By Vector
// **********************************************************************
            
      interface rsbv_dia
        module procedure irsbv_dia
        module procedure srsbv_dia
        module procedure drsbv_dia
        module procedure crsbv_dia
        module procedure zrsbv_dia
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void irsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,j,n,lda,ndiag;
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=0;
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                    x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=0;
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
            }else{/*else*/ if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  irsbv_dia 
// **********************************************************************
// **********************************************************************
      void srsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,j,n,lda,ndiag;
      char  diag,part;
      float de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=0.0e0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0.0e0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                    x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=0.0e0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
            }else{/*else*/ if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0.0e0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  srsbv_dia 
// **********************************************************************
// **********************************************************************
      void drsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,j,n,lda,ndiag;
      char  diag,part;
     double de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=0.0d0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0.0d0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                    x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=0.0d0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
            }else{/*else*/ if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=0.0d0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  drsbv_dia 
// **********************************************************************
// **********************************************************************
      void crsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,lda,ndiag;
      char  diag,part;
      comple_f de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=(0.0e0,0.0e0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=(0.0e0,0.0e0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                    x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=(0.0e0,0.0e0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
            }else{/*else*/ if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=(0.0e0,0.0e0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  crsbv_dia 
// **********************************************************************
// **********************************************************************
      void zrsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,lda,ndiag;
      char  diag,part;
      complex_d de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=(0.0d0,0.0d0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((mat->IA1[j]>-i)&&(mat->IA1[j]<0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=(0.0d0,0.0d0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                    x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=(0.0d0,0.0d0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
            }else{/*else*/ if((mat->IA1[j]<n-i+1)&&(mat->IA1[j]>0)) {//then
                     x[i]-=mat->A[lda*j+i]*x[i+mat->IA1[j]];
                  }//end if
               }//end do
               if (de!=(0.0d0,0.0d0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zrsbv_dia 
// **********************************************************************
// **********************************************************************
      end module mod_rsbv_dia
