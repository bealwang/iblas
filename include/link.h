#ifndef spblas_c_link_h
#define spblas_c_link_h
#include <malloc.h>
#include "comm_tools.h"
#include "types.h"
#include "properties.h"
// **********************************************************************
//     Author : luoyulong
//     Date of last modification : 8.29.2012
//     Description : THE PRINCIPAL DATA STRUCTURE
//
//                   new: creates new node WITHOUT initialization
//                   del: frees unused memory, does NOT care if there
//                        is other memory that should be freed first
//                   accessdata:  returns a pointer to the matrix
//                        inside the relevant node
// **********************************************************************

typedef struct isp_l
{
  ISPMAT contents;
  int number;
  struct isp_l* pntr;
}isp_linknode;

typedef struct ssp_l
{
  SSPMAT contents;
  int number;
 struct ssp_l* pntr;
}ssp_linknode;

typedef struct dsp_l
{
  DSPMAT contents;
  int number;
  struct dsp_l* pntr;
}dsp_linknode;

typedef struct csp_l
{
  CSPMAT contents;
  int number;
 struct csp_l* pntr;
}csp_linknode;

typedef struct zsp_l
{
  ZSPMAT contents;
  int number;
  struct zsp_l* pntr;
}zsp_linknode;

isp_linknode* new_isp(int *ierr);

void del_isp(isp_linknode* sp_l,int *ierr);

ISPMAT * accessdata_isp(isp_linknode* sp_l,int *ierr);

ssp_linknode * new_ssp(int *ierr);

void del_ssp(ssp_linknode * sp_l,int *ierr);

SSPMAT * accessdata_ssp(ssp_linknode* sp_l,int *ierr);

dsp_linknode * new_dsp(int *ierr);

void del_dsp(dsp_linknode * sp_l,int *ierr);

DSPMAT * accessdata_dsp(dsp_linknode* sp_l,int *ierr);

csp_linknode * new_csp(int *ierr);

void del_csp(csp_linknode * sp_l,int *ierr);

CSPMAT * accessdata_csp(csp_linknode* sp_l,int *ierr);

zsp_linknode * new_zsp(int *ierr);

void del_zsp(zsp_linknode * sp_l,int *ierr);

ZSPMAT * accessdata_zsp(zsp_linknode* sp_l,int *ierr);
#endif
