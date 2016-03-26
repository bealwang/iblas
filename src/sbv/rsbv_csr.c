      module mod_rsbv_csr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'CSR'-STORAGE
//                   rsbv=Right Solve By Vector
// **********************************************************************
            
      interface rsbv_csr
        module procedure irsbv_csr
        module procedure srsbv_csr
        module procedure drsbv_csr
        module procedure crsbv_csr
        module procedure zrsbv_csr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void irsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,base,ofs,pntr
      char  diag,part;
      int de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0; 
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  irsbv_csr 
// **********************************************************************
// **********************************************************************
      void srsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,base,ofs,pntr
      char  diag,part;
      float de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0.0e0  
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0.0e0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0.0e0 
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0.0e0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  srsbv_csr 
// **********************************************************************
// **********************************************************************
      void drsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,base,ofs,pntr
      char  diag,part;
     double de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0.0d0  
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0.0d0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=0.0d0 
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==0.0d0 ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  drsbv_csr 
// **********************************************************************
// **********************************************************************
      void crsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,pntr
      char  diag,part;
      comple_f de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=(0.0e0,0.0e0)  
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==(0.0e0,0.0e0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=(0.0e0,0.0e0) 
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==(0.0e0,0.0e0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  crsbv_csr 
// **********************************************************************
// **********************************************************************
      void zrsbv_csr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,ofs,pntr
      char  diag,part;
      complex_d de;
      *ierr=-1;
      n=n_x;
      if ((mat->FIDA!='CSR')||(mat->M!=n)||(mat->K!=n)) {//then
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
         for(i=0;i<n;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=(0.0d0,0.0d0)  
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==(0.0d0,0.0d0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }else{/*else*/ 
         for(i=n-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
               pntr++;
               }//end do
            }else{/*else*/
               pntr=mat->pb[i];
               de=(0.0d0,0.0d0) 
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     x[i]=x[i]                  -mat->A[pntr+ofs]*x(mat->IA1[pntr+ofs]+ofs) 
                  }else{/*else*/
                     de=mat->A[pntr+ofs];
                  }//end if
                  pntr++;
               }//end do
               if(de==(0.0d0,0.0d0) ) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  x[i]=x[i]/de;
               }//end if
            }//end if
         }//end do
         *ierr=0;
      }//end if
      }//function end  zrsbv_csr 
// **********************************************************************
// **********************************************************************
      end module mod_rsbv_csr
