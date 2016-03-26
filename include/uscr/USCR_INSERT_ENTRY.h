#ifndef USCR_INSERT_ENTRY_H
#define USCR_INSERT_ENTRY_H
#include "blas_enum.h"
#include "INSERTING.h"
#include "blas_sparse_namedconstants.h"
#include "INS_ROUTINER.h"
//int BLAS_iuscr_insert_entry(i_matrix* A,SCALAR_IN val,int i,int j);
//int BLAS_suscr_insert_entry(blas_sparse_matrix A,SCALAR_IN val,int i,int j);
int BLAS_duscr_insert_entry(d_matrix* A,double val,int i,int j);
//int BLAS_cuscr_insert_entry(blas_sparse_matrix A,SCALAR_IN val,int i,int j);
//int BLAS_zuscr_insert_entry(blas_sparse_matrix A,SCALAR_IN val,int i,int j);


//// **********************************************************************
//int BLAS_iuscr_insert_entry(i_matrix* A,int val,int i,int j)
//  {
//  d_matrix* pmatrix;

//      pmatrix=daccess_matrix(A);

//      switch (pmatrix->format)
//      {
//       case 'b':
//          iINS_bl_entr (pmatrix,val,i,j);
//          break;
//      case 'v':
//           iINS_varbl_entr (pmatrix,val,i,j,istat);
//           break;
//      case 'n':
//          iINS_entry (pmatrix,val,i,j,istat);
//          break;
//      default:
//         istat = blas_error_param;;
//         return;
//      }
//      }//end subroutine iuscr_insert_entry
// **********************************************************************
// **********************************************************************

//int  BLAS_suscr_insert_entry  (s_matrix* A,int val,int i,int j)
//      {
//        d_matrix* pmatrix;
//        istat=-1;
//         pmatrix=daccess_matrix(A);
//        if (istat!=0) return;
//        switch (pmatrix->format)
//        {
//          case 'b':
//            dINS_bl_entr (pmatrix,val,i,j,istat);
//          case'v':
//             dINS_varbl_entr(pmatrix,val,i,j,istat);
//          case 'n':
//            dINS_entry (pmatrix,val,i,j,istat);
//          default:
//           istat = blas_error_param;
//          return;
//        }
//      }//end subroutine suscr_insert_entry
//// **********************************************************************
//////***************************************************************************
//// **********************************************************************
//      int BLAS_duscr_insert_entry  (d_matrix* A,double val,int i,int j)
//      {
//      d_matrix* pmatrix;
//      int istat=-1;
//       pmatrix=daccess_matrix(A);
//      switch (pmatrix->format)
//      {
//        case 'b':
//          dINS_bl_entr (pmatrix,val,i,j,&istat);
//        case'v':
//           dINS_varbl_entr(pmatrix,val,i,j,&istat);
//        case 'n':
//          dINS_entry (pmatrix,val,i,j,&istat);
//        default:
//         istat=blas_error_param;
//        return;
//      }
//      return istat;
//      }//end subroutine duscr_insert_entry
// **********************************************************************
////***************************************************************************
// **********************************************************************
//      subroutine cuscr_insert_entry  (a,val,i,j,istat)
//      {
//      complex(KIND=sp)  ,intent(in) ::val
//      integer ,intent(in) ::a,i,j
//      integer,intent(out)::istat
//      type(c_matrix ),pointer ::pmatrix
//      istat=-1
//       access_matrix(pmatrix,a,istat)
//      if (istat!=0) return
//      switch (pmatrix->format)
//      case('b')
//          cINS_bl_entr (pmatrix,val,i,j,istat)
//      case('v')
//           cINS_varbl_entr (pmatrix,val,i,j,istat)
//      case('n')
//          cINS_entry (pmatrix,val,i,j,istat)
//      case default
//         istat = blas_error_param
//         return
//      end select
//      }//end subroutine cuscr_insert_entry
//// **********************************************************************
//////***************************************************************************
//// **********************************************************************
//      subroutine zuscr_insert_entry  (a,val,i,j,istat)
//      {
//      complex(KIND=dp)  ,intent(in) ::val
//      integer ,intent(in) ::a,i,j
//      integer,intent(out)::istat
//      type(z_matrix ),pointer ::pmatrix
//      istat=-1
//       daccess_matrix(pmatrix,a,istat)
//      if (istat!=0) return
//      switch (pmatrix->format)
//      case('b')
//          zINS_bl_entr (pmatrix,val,i,j,istat)
//      case('v')
//           zINS_varbl_entr (pmatrix,val,i,j,istat)
//      case('n')
//          zINS_entry (pmatrix,val,i,j,istat)
//      case default
//         istat = blas_error_param
//         return
//      end select
//      }//end subroutine zuscr_insert_entry
//// **********************************************************************
//////***************************************************************************
//// **********************************************************************
//end module mod_uscr_insert_entry





#endif // USCR_INSERT_ENTRY_H
