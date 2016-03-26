#ifndef USCONV_H
#define USCONV_H
void usconv_bco2bdi(void* spmat,int* ierr);
void usconv_bco2bsc(void* spmat,int* ierr);
void usconv_bco2bsr(void* spmat,int* ierr);
void usconv_bco2coo(void* spmat,int* ierr);
void usconv_bdi2bco(void* spmat,int* ierr);
void usconv_bsc2bco(void* spmat,int* ierr);
void usconv_bsr2bco(void* spmat,int* ierr);
void usconv_coo2csr(void* spmat,int* ierr);
void usconv_coo2csc(void* spmat,int* ierr);
void usconv_coo2dia(void* spmat,int* ierr);
void usconv_coo2bco(void* spmat,int* ierr);
void usconv_csc2coo(void* spmat,int* ierr);
void usconv_csr2coo(void* spmat,int* ierr);
void usconv_csr2bco(void* spmat,int* ierr);
void usconv_dia2coo(void* spmat,int* ierr);
#endif // USCONV_H
