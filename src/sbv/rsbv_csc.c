      module mod_rsbv_csc
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSC'-STORAGE
//                   rsbv=Right Solve By Vector
// **********************************************************************
            
      interface rsbv_csc
        module procedure irsbv_csc
        module procedure srsbv_csc
        module procedure drsbv_csc
        module procedure crsbv_csc
        module procedure zrsbv_csc
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void irsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  j,n,base,ofs,pntr;
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  irsbv_csc 
// **********************************************************************
// **********************************************************************
      void srsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  j,n,base,ofs,pntr;
      char  diag,part;
      float de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0.0e0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0.0e0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  srsbv_csc 
// **********************************************************************
// **********************************************************************
      void drsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  j,n,base,ofs,pntr;
      char  diag,part;
     double de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0.0d0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==0.0d0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  drsbv_csc 
// **********************************************************************
// **********************************************************************
      void crsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,ofs,pntr;
      char  diag,part;
      comple_f de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==(0.0e0,0.0e0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==(0.0e0,0.0e0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  crsbv_csc 
// **********************************************************************
// **********************************************************************
      void zrsbv_csc (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,ofs,pntr;
      char  diag,part;
      complex_d de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSC')||(mat->M!=n)||(mat->K!=n)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
       get_descra(mat->DESCRA,'d',&diag,*ierr);
      
       get_descra(mat->DESCRA,'a',&part,*ierr);
      
      if ((part!='U')&&(part!='L')) {//then
         *ierr=blas_error_param;
         return;
      }//end if
      *ierr=-1;
      if (part=='L') {//then
         for(j=0;j<n;j++){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==(0.0d0,0.0d0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/
         for(j=n-1;j>=0;j--){//for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*x[j] ; 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[j];
               while((pntr<mat->pe[j])&&((mat->IA1[pntr+ofs]+ofs)!=j)){//while
                  pntr++;
               }//end do
               if((mat->IA1[pntr+ofs]+ofs)==j) {//then
                  de=mat->A[pntr+ofs];
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(de==(0.0d0,0.0d0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[j]=x[j]/de;
                  de=x[j];
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  x[mat->IA1[pntr+ofs]+ofs]-=mat->A[pntr+ofs]*de;
                  pntr++;
               }//end do
               x[j]=de;
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zrsbv_csc 
// **********************************************************************
// **********************************************************************
      end module mod_rsbv_csc
