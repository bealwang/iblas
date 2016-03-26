#ifndef BLAS_ENUM_H
#define BLAS_ENUM_H

typedef enum trans_type
{
  ORIGIN_MATRIX,
  TRANSP_MATRIX
}blas_trans_type;
typedef enum store_format
{
  COO_FORMAT,
  CSR_FORMAT,
  CSC_FORMAT,
  DIA_FORMAT,
  BCO_FORMAT,
  BSR_FORMAT,
  BSC_FORMAT,
  BDI_FORMAT,
  VBR_FORMAT,
  NULL_FORMAT
}blas_store_format;
enum blas_base_type
{
  BLAS_BASE_ZERO=0,
  BLAS_BASE_ONE=1
};
typedef enum order_type
{
  type_01,
  type_02
}blas_order_type;
typedef struct
{

} blas_sparse_matrix ;


#endif // BLAS_ENUM_H
