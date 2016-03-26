      module mod_lsbv_csr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH TRANSPOSE IN 'CSR'-STORAGE
//                   lsbv = Left Solve By Vector
// **********************************************************************      
      use representation_of_data
      use properties
     
      interface lsbv_csr
        module procedure ilsbv_csr
        module procedure slsbv_csr
        module procedure dlsbv_csr
        module procedure clsbv_csr
        module procedure zlsbv_csr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  j,n,base,ofs,pntr;
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  ilsbv_csr 
// **********************************************************************
// **********************************************************************
      void slsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  j,n,base,ofs,pntr;
      char  diag,part;
      float de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0.0e0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0.0e0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  slsbv_csr 
// **********************************************************************
// **********************************************************************
      void dlsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  j,n,base,ofs,pntr;
      char  diag,part;
     double de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0.0d0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== 0.0d0 ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  dlsbv_csr 
// **********************************************************************
// **********************************************************************
      void clsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,ofs,pntr;
      char  diag,part;
      comple_f de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== (0.0e0, 0.0e0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== (0.0e0, 0.0e0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  clsbv_csr 
// **********************************************************************
// **********************************************************************
      void zlsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,ofs,pntr;
      char  diag,part;
      complex_d de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== (0.0d0, 0.0d0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                 x[mat->IA1[pntr + ofs] + ofs]-=mat->A[pntr + ofs] * x[j]; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&& (mat->IA1[pntr + ofs] + ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr + ofs] + ofs==j) {//then
                  de=mat->A[pntr + ofs];
               }else{/*else*/
                  *ierr = blas_error_singtria;
                  return;
               }//end if
               if(de== (0.0d0, 0.0d0) ) {//then
                  *ierr = blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j] = x[j]/de;
                  de = x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr + ofs] + ofs]-= mat->A[pntr + ofs] * de;
                  pntr++;
               }//end do
               x[j] = de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zlsbv_csr 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_csr
