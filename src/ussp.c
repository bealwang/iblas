/*
 * m is property, 
 * a is number 
 */
// a doesnot change at all ,what is the fuction of this ?
int ussp(INOUT int *a, IN int m, OUT int *istat)
{
    int b, rest;
    i_matrix * ipmatrix;
    d_matrix * dpmatrix;
    s_matrix * spmatrix;
    c_matrix * cpmatrix;
    z_matrix * zpmatrix;
    b = -a; //a is negetive number
    *istat = 0;
    rest =  b%no_of_types;
    switch (rest)
    {
        case(ISP_MATRIX)
            access_matrix(ipmtrix,a, istat);
            ipmatrix->property = m; // what is m mean ?
            break;
        case(SSP_MATRIX)
            access_matrix(spmtrix,a, istat);
            spmatrix->property = m; // what is m mean ?
            break;
        case(DSP_MATRIX)
            access_matrix(dpmtrix,a, istat);
            dpmatrix->property = m; // what is m mean ?
            break;
        case(CSP_MATRIX)
            access_matrix(cpmtrix,a, istat);
            cpmatrix->property = m; // what is m mean ?
            break;
        case(ZSP_MATRIX)
            access_matrix(zpmtrix,a, istat);
            zpmatrix->property = m; // what is m mean ?
            break;
        case default:
            *istat = blas_error_para
            return -1;
    }
    
}
