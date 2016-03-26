      module mod_rsbv_vbr
// **********************************************************************
//     Author : C. Voemel
//     Date of last modification : 7.7.00
//     Description : PERFORMS TRI. SOLVE WITH MATRIX IN 'VBR'-STORAGE
//                   rsbv=Right Solve By Vector
// **********************************************************************
                  
      interface rsbv_vbr
        module procedure irsbv_vbr
        module procedure srsbv_vbr
        module procedure drsbv_vbr
        module procedure crsbv_vbr
        module procedure zrsbv_vbr
      end interface
      #include "sbv.h"
// **********************************************************************
// **********************************************************************
      void irsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      int  i,n,base,pntr,ofs,mb,nb,dd;
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
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function irsbv_vbr 
// **********************************************************************
// **********************************************************************
      void srsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
            int  i,n,base,pntr,ofs,mb,nb,dd;
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
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function srsbv_vbr 
// **********************************************************************
// **********************************************************************
      void drsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
      
      
            int  i,n,base,pntr,ofs,mb,nb,dd;
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
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function drsbv_vbr 
// **********************************************************************
// **********************************************************************
      void crsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,pntr,ofs,mb,nb,dd;
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
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function crsbv_vbr 
// **********************************************************************
// **********************************************************************
      void zrsbv_vbr (ISPMAT* mat,int* x,int n_x,int* ierr)
      {
            int  i,n,base,pntr,ofs,mb,nb,dd;
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
         for(i=0;i<mb;i++){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_left_lower(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }else{/*else*/ 
         for(i=mb-1;i>=0;i--){//for
            if (diag=='U') {//then
               pntr=mat->pb[i];
               while(pntr<mat->pe[i]){//while
                  start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                  end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_y=mat->bp1[i]+ofs;
                  end_y=mat->bp1[i+1]+ofs-1;
                  len_y=end_y-start_y+1;
                  start_a=mat->IA2[pntr+ofs]+ofs;
                  end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                   block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  pntr++;
               }//end do 
            }else{/*else*/
               pntr=mat->pb[i];
               dd=-1;
               while(pntr<mat->pe[i]){//while
if((mat->IA1[pntr+ofs]+ofs)!=i) {//then

                     start_x=mat->bp2[mat->IA1[pntr+ofs]+ofs]+ofs;
                     end_x=mat->bp2[mat->IA1[pntr+ofs]+ofs+1]+ofs-1;
                     len_x=end_x-start_x+1;
                     start_y=mat->bp1[i]+ofs;
                     end_y=mat->bp1[i+1]+ofs-1;
                     len_y=end_y-start_y+1;
                     start_a=mat->IA2[pntr+ofs]+ofs;
                     end_a=mat->IA2[pntr+ofs+1]+ofs-1;
                      block_mult_vec(mat->A[start_a][end_a],x[start_x][end_x],len_x,y[start_y][end_y],len_y,store,*ierr); 
             x[start_y][end_y]-=y[start_y][end_y];
                  }else{/*else*/
                     dd=pntr;
                  }//end if 
                  pntr++;
               }//end do 
               if(dd==-1) {//then
                  *ierr=blas_error_singtria;
                  return;
               }else{/*else*/
                  start_x=mat->bp1[i]+ofs;
                  end_x=mat->bp1[i+1]+ofs-1;
                  len_x=end_x-start_x+1;
                  start_a=mat->IA2[dd+ofs]+ofs;
                  end_a=mat->IA2[dd+ofs+1]+ofs-1;
                   invert_right_upper(mat->A[start_a][end_a],x[start_x][end_x],len_x,store,*ierr);
               }//end if
               if(*ierr!=0) {//then
                  *ierr=blas_error_singtria;
                  return;
               }//end if
            }//end if
         }//end do
         free(y);
         if(*ierr!=0) {//then
            *ierr=blas_error_memdeloc;
            return;
         }//end if
      }//end if
} //end function zrsbv_vbr 
// **********************************************************************
// **********************************************************************
      end module mod_rsbv_vbr
