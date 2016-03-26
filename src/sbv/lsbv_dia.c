      module mod_lsbv_dia
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'DIA'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************
                 
      interface lsbv_dia
        module procedure ilsbv_dia
        module procedure slsbv_dia
        module procedure dlsbv_dia
        module procedure clsbv_dia
        module procedure zlsbv_dia
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,j,n,lda,ndiag;
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0)) {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de=0;
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0))                   {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de=0;
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0))                   {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  ilsbv_dia 
// **********************************************************************
// **********************************************************************
      void slsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,j,n,lda,ndiag;
      char  diag,part;
      float de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0)) {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de = 0.0e0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0))                   {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0.0e0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de = 0.0e0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0))                   {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0.0e0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  slsbv_dia 
// **********************************************************************
// **********************************************************************
      void dlsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,j,n,lda,ndiag;
      char  diag,part;
     double de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0)) {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de = 0.0d0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0))                   {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0.0d0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de = 0.0d0 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0))                   {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= 0.0d0 ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  dlsbv_dia 
// **********************************************************************
// **********************************************************************
      void clsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,lda,ndiag;
      char  diag,part;
      comple_f de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0)) {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de = (0.0e0,0.0e0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0))                   {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= (0.0e0,0.0e0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de = (0.0e0,0.0e0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0))                   {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= (0.0e0,0.0e0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  clsbv_dia 
// **********************************************************************
// **********************************************************************
      void zlsbv_dia (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,j,n,lda,ndiag;
      char  diag,part;
      complex_d de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='DIA')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr = blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'d',&lda,*ierr);
      
       get_infoa(mat->INFOA,'e',&ndiag,*ierr);
      
      *ierr=-1;
      if (part=='L') {//then
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                 if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0)) {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                 }//end if
               }//end do
            }else{/*else*/
               de = (0.0d0,0.0d0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
               }else{/*else*/ if((-mat->IA1[j]<n-i+1)&&(mat->IA1[j]<0))                   {//then
             x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= (0.0d0,0.0d0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               for(j=0;j<n;j++){//fordiag
                  if((mat->IA1[j]<i)&&(mat->IA1[j]>0)) {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
            }else{/*else*/
               de = (0.0d0,0.0d0) 
               for(j=0;j<n;j++){//fordiag
                  if (mat->IA1[j]==0) {//then
                     de=mat->A[lda*j+i];
                  }else{/*else*/ if((mat->IA1[j]<i)&&(mat->IA1[j]>0))                   {//then
            x[i]-=mat->A[lda*j+i-mat->IA1[j]]*x[i-mat->IA1[j]];
                  }//end if
               }//end do
               if (de!= (0.0d0,0.0d0) ) {//then 
                  x[i]=x[i]/de;
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zlsbv_dia 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_dia
