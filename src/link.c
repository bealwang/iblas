#include "link.h"
isp_linknode* new_isp(int *ierr)
{
  *ierr=-1;
  isp_linknode * sp_l=(isp_linknode *)aligned_malloc (sizeof(isp_linknode));
  sp_l->number=ISP_MATRIX;

  sp_l->contents.A=NULL;
  sp_l->contents.IA1=NULL;
  sp_l->contents.IA2=NULL;
  sp_l->contents.PB=NULL;
  sp_l->contents.PE=NULL;
  sp_l->contents.BP1=NULL;
  sp_l->contents.BP2=NULL;

  //sp_l->contents.FIDA=' ';
  //sp_l->contents.DESCRA=' ';
  //sp_l->contents.INFOA=0;

  sp_l->pntr=NULL;
return sp_l;
}
void del_isp(isp_linknode* sp_l,int *ierr)
{
  *ierr=0;
  if(sp_l!=NULL) aligned_free (sp_l);
  else *ierr=-1;
}

ISPMAT * accessdata_isp(isp_linknode* sp_l,int *ierr)
{
  if(sp_l!=NULL)
    {
      *ierr=0;
      return &(sp_l->contents);
    }
  else
    {
      *ierr=-1;
      return NULL;
    }
}



ssp_linknode * new_ssp(int *ierr)
{
  *ierr=-1;
  ssp_linknode * sp_l=(ssp_linknode *)aligned_malloc(sizeof(ssp_linknode));
  sp_l->number=SSP_MATRIX;

  sp_l->contents.A=NULL;
  sp_l->contents.IA1=NULL;
  sp_l->contents.IA2=NULL;
  sp_l->contents.PB=NULL;
  sp_l->contents.PE=NULL;
  sp_l->contents.BP1=NULL;
  sp_l->contents.BP2=NULL;

 // sp_l->contents.FIDA='';
 // sp_l->contents.DESCRA='';
  //sp_l->contents.INFOA=0;

  sp_l->pntr=NULL;
return sp_l;
}
void del_ssp(ssp_linknode * sp_l,int *ierr)
{
  *ierr=0;
  if(sp_l!=NULL) aligned_free(sp_l);
  else *ierr=-1;
}

SSPMAT * accessdata_ssp(ssp_linknode* sp_l,int *ierr)
{
  if(sp_l!=NULL)
    {
      *ierr=0;
      return &(sp_l->contents);
    }
  else
    {
      *ierr=-1;
      return NULL;
    }
}

dsp_linknode * new_dsp(int *ierr)
{
  int i=0;
  *ierr=-1;
  dsp_linknode * sp_l=(dsp_linknode *)aligned_malloc(sizeof(dsp_linknode));
  sp_l->number=DSP_MATRIX;

  sp_l->contents.A=NULL;
  sp_l->contents.IA1=NULL;
  sp_l->contents.IA2=NULL;
  sp_l->contents.PB=NULL;
  sp_l->contents.PE=NULL;
  sp_l->contents.BP1=NULL;
  sp_l->contents.BP2=NULL;

  sp_l->contents.FIDA=NULL_FORMAT;

  for(i=0;i<11;i++)
    sp_l->contents.DESCRA[i]='0';

  for(i=0;i<10;i++)
  sp_l->contents.INFOA[i]='0';

  sp_l->pntr=NULL;
  *ierr=0;
return sp_l;
}
void del_dsp(dsp_linknode * sp_l,int *ierr)
{
  *ierr=0;
  if(sp_l!=NULL) aligned_free(sp_l);
  else *ierr=-1;
}

DSPMAT * accessdata_dsp(dsp_linknode* sp_l,int *ierr)
{
  if(sp_l!=NULL)
    {
      *ierr=0;
      return &(sp_l->contents);
    }
  else
    {
      *ierr=-1;
      return NULL;
    }
}

csp_linknode * new_csp(int *ierr)
{
  *ierr=-1;
  csp_linknode * sp_l=(csp_linknode *)aligned_malloc (sizeof(csp_linknode));
  sp_l->number=CSP_MATRIX;

  sp_l->contents.A=NULL;
  sp_l->contents.IA1=NULL;
  sp_l->contents.IA2=NULL;
  sp_l->contents.PB=NULL;
  sp_l->contents.PE=NULL;
  sp_l->contents.BP1=NULL;
  sp_l->contents.BP2=NULL;

 // sp_l->contents.FIDA='';
  //sp_l->contents.DESCRA='';
  //sp_l->contents.INFOA=0;

  sp_l->pntr=NULL;
return sp_l;
}
void del_csp(csp_linknode * sp_l,int *ierr)
{
  *ierr=0;
  if(sp_l!=NULL) aligned_free(sp_l);
  else *ierr=-1;
}

CSPMAT * accessdata_csp(csp_linknode* sp_l,int *ierr)
{
  if(sp_l!=NULL)
    {
      *ierr=0;
      return &(sp_l->contents);
    }
  else
    {
      *ierr=-1;
      return NULL;
    }
}

zsp_linknode * new_zsp(int *ierr)
{
  *ierr=-1;
  zsp_linknode * sp_l=(zsp_linknode *)aligned_malloc (sizeof(zsp_linknode));
  sp_l->number=ZSP_MATRIX;

  sp_l->contents.A=NULL;
  sp_l->contents.IA1=NULL;
  sp_l->contents.IA2=NULL;
  sp_l->contents.PB=NULL;
  sp_l->contents.PE=NULL;
  sp_l->contents.BP1=NULL;
  sp_l->contents.BP2=NULL;

 // sp_l->contents.FIDA='';
 // sp_l->contents.DESCRA='';
  //sp_l->contents.INFOA=0;

  sp_l->pntr=NULL;
return sp_l;
}
void del_zsp(zsp_linknode * sp_l,int *ierr)
{
  *ierr=0;
  if(sp_l!=NULL) aligned_free(sp_l);
  else *ierr=-1;
}

ZSPMAT * accessdata_zsp(zsp_linknode* sp_l,int *ierr)
{
  if(sp_l!=NULL)
    {
      *ierr=0;
      return &(sp_l->contents);
    }
  else
    {
      *ierr=-1;
      return NULL;
    }
}
