      module mod_lsbv_vbr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'VBR'-STORAGE
//                   lsbv=Left Solve By Vector
// **********************************************************************
                  
      interface lsbv_vbr
        module procedure ilsbv_vbr      
        module procedure slsbv_vbr
        module procedure dlsbv_vbr
        module procedure clsbv_vbr
        module procedure zlsbv_vbr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void ilsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  j,n,base,pntr,ofs,mb,nb,dd;
      int  start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y;
      char diag,part,store;
      int* y;
      *ierr=-1;
      n=n_x;
      y=malloc(sizeof(y*)*n_x);
      
      y=0 
      if ((mat->FIDA!='VBR')||(mat->M!=n)||(mat->K!=n)||(size(mat->bp1)!=size(mat->bp2))||(maxval(abs(mat->bp1-mat->bp2))!=0)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
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
      start_a=-1;
      end_a=-1;
      start_x=-1;
      end_x=-1;
      start_y=-1;
      end_y=-1;
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function ilsbv_vbr 
// **********************************************************************
// **********************************************************************
      void slsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  j,n,base,pntr,ofs,mb,nb,dd;
      int  start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y;
      char diag,part,store;
      float* y;
      *ierr=-1;
      n=n_x;
      y=malloc(sizeof(y*)*n_x);
      
      y=0.0e0 
      if ((mat->FIDA!='VBR')||(mat->M!=n)||(mat->K!=n)||(size(mat->bp1)!=size(mat->bp2))||(maxval(abs(mat->bp1-mat->bp2))!=0)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
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
      start_a=-1;
      end_a=-1;
      start_x=-1;
      end_x=-1;
      start_y=-1;
      end_y=-1;
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function slsbv_vbr 
// **********************************************************************
// **********************************************************************
      void dlsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  j,n,base,pntr,ofs,mb,nb,dd;
      int  start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y;
      char diag,part,store;
      double* y;
      *ierr=-1;
      n=n_x;
      y=malloc(sizeof(y*)*n_x);
      
      y=0.0d0 
      if ((mat->FIDA!='VBR')||(mat->M!=n)||(mat->K!=n)||(size(mat->bp1)!=size(mat->bp2))||(maxval(abs(mat->bp1-mat->bp2))!=0)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
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
      start_a=-1;
      end_a=-1;
      start_x=-1;
      end_x=-1;
      start_y=-1;
      end_y=-1;
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function dlsbv_vbr 
// **********************************************************************
// **********************************************************************
      void clsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,pntr,ofs,mb,nb,dd;
      int  start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y;
      char diag,part,store;
      complex_f* y;
      *ierr=-1;
      n=n_x;
      y=malloc(sizeof(y*)*n_x);
      
      y=(0.0e0,0.0e0) 
      if ((mat->FIDA!='VBR')||(mat->M!=n)||(mat->K!=n)||(size(mat->bp1)!=size(mat->bp2))||(maxval(abs(mat->bp1-mat->bp2))!=0)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
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
      start_a=-1;
      end_a=-1;
      start_x=-1;
      end_x=-1;
      start_y=-1;
      end_y=-1;
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function clsbv_vbr 
// **********************************************************************
// **********************************************************************
      void zlsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  j,n,base,pntr,ofs,mb,nb,dd;
      int  start_a,end_a,start_x,end_x,len_x,start_y,end_y,len_y;
      char diag,part,store;
      complex_d* y;
      *ierr=-1;
      n=n_x;
      y=malloc(sizeof(y*)*n_x);
      
      y=(0.0d0,0.0d0) 
      if ((mat->FIDA!='VBR')||(mat->M!=n)||(mat->K!=n)||(size(mat->bp1)!=size(mat->bp2))||(maxval(abs(mat->bp1-mat->bp2))!=0)) {//then
         *ierr=blas_error_param;
         return;
      }//end if
       get_infoa(mat->INFOA,'b',&base,*ierr);
      
      ofs=1-base;
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
      start_a=-1;
      end_a=-1;
      start_x=-1;
      end_x=-1;
      start_y=-1;
      end_y=-1;
      if (part=='L') {//then
         for(j=mb-1;j>=0;j--){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(j=0;j<mb;j++){//end for
            if (diag=='U') {//then
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[j];
               dd=-1;
               while((pntr<mat->pe[j])&&(mat->IA1[pntr+ofs]+ofs!=j)){//while
                  pntr++;
               }//end do
               if(mat->IA1[pntr+ofs]+ofs==j) {//then
                  dd=pntr;
               }else{/*else*/
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[j]+ofs;
                  end_x=mat->bp1[j+1]+ofs -1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   invert_T_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
               pntr=mat->pb[j];
               while(pntr<mat->pe[j]){//while
                  if(mat->IA1[pntr+ofs]+ofs!=j) {//then
                     start_x=mat->bp1[j]+ofs;
                     end_x=mat->bp1[j+1]+ofs -1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_y=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_T_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
               x[start_y][end_y]-=y[start_y][end_y];
                  }//end if
                  pntr++;
               }//end do 
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function zlsbv_vbr 
// **********************************************************************
// **********************************************************************
      end module mod_lsbv_vbr
