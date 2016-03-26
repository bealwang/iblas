#ifndef INS_ROUTINER_H
#define INS_ROUTINER_H
#include "INSERTING.h"
#include "comm_tools.h"
#include "uscr/uscr_bco.h"
#include "uscr/uscr_vbr.h"
#include "uscr/uscr_coo.h"
void dINS_entry(d_matrix * pmatrix, double val,int i,int j,int* istat);

void dINS_varbl_entr(d_matrix* vpmatrix,double val,int i,int j,int* istat);

void dINS_varblock(d_matrix* vpmatrix,double** val,int r_val,int c_val,int i,int j,int* istat);

void dINS_block (d_matrix* pmatrix,double** val,int r_val,int c_val,int i,int j,int* istat);

void dINS_bl_entr (d_matrix* pmatrix,double val,int i,int j,int* istat);

dsp_linknode* duscr_blockend(d_matrix* A,int prpty,int *istat);

dsp_linknode*  duscr_varend (d_matrix* A,int prpty,int *istat);

dsp_linknode* duscr_normend(d_matrix* A,int prpty,int *istat);
#endif
