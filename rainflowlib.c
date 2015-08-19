#include <math.h>

/* ++++++++++ BEGIN RF3 [ampl ampl_mean nr_of_cycle] */
/* ++++++++++ Rain flow without time analysis */
//Copyright (c) 1999-2002 by Adam Nieslony
//Visit the MATLAB Central File Exchange for latest version
//http://www.mathworks.com/matlabcentral/fileexchange/3026
void rf3(double *array_ext, double *array_out, int n, int *nout) {
    
    double *pr, *po, a[16384], ampl, mean;
    int tot_num, index, j, cNr1, cNr2;    
    
    tot_num = n;
    
    // pointers to the first element of the arrays
    pr = &array_ext[0];
    po = &array_out[0];
    
    // The original rainflow counting by Nieslony, unchanged
    j = -1;
    cNr1 = 1;
    for (index=0; index<tot_num; index++) {
        a[++j]=*pr++;
        while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
            ampl=fabs( (a[j-1]-a[j-2])/2 );
            switch(j) {
                case 0: { break; }
                case 1: { break; }
                case 2: {
                    mean=(a[0]+a[1])/2;
                    a[0]=a[1];
                    a[1]=a[2];
                    j=1;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=0.50;
                    }
                    break;
                }
                default: {
                    mean=(a[j-1]+a[j-2])/2;
                    a[j-2]=a[j];
                    j=j-2;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=1.00;
                        cNr1++;
                    }
                    break;
                }
            }
        }
    }
    cNr2 = 1;
    for (index=0; index<j; index++) {
        ampl=fabs(a[index]-a[index+1])/2;
        mean=(a[index]+a[index+1])/2;
        if (ampl > 0){
            *po++=ampl;
            *po++=mean;
            *po++=0.50;
            cNr2++;
        }
    }
    // array of ints nout is outputted
    nout[0] = cNr1;
    nout[1] = cNr2;
}
/* ++++++++++ END RF3 */


// ++ BEGIN RF5 [ampl ampl_mean nr_of_cycle cycle_begin_time cycle_period_time]
/* ++++++++++ Rain flow with time analysis */
//Copyright (c) 1999-2002 by Adam Nieslony
//Visit the MATLAB Central File Exchange for latest version
//http://www.mathworks.com/matlabcentral/fileexchange/3026
void
rf5(double *array_ext, double *array_t, double *array_out, int n, int *nout) {
    double *pr, *pt, *po, a[16384], t[16384], ampl, mean, period, atime;
    int tot_num, index, j, cNr1, cNr2;
    
    
//    tot_num = mxGetM(array_ext) * mxGetN(array_ext);
    tot_num = n;

    // pointers to the first element of the arrays
    pr = &array_ext[0];
    pt = &array_t[0];
    po = &array_out[0];
    
//    array_out = mxCreateDoubleMatrix(5, tot_num-1, mxREAL);
    
    // The original rainflow counting by Nieslony, unchanged
    j = -1;
    cNr1 = 1;
    for (index=0; index<tot_num; index++) {
        a[++j]=*pr++;
        t[j]=*pt++;
        while ( (j >= 2) && (fabs(a[j-1]-a[j-2]) <= fabs(a[j]-a[j-1])) ) {
            ampl=fabs( (a[j-1]-a[j-2])/2 );
            switch(j)
{
                case 0: { break; }
                case 1: { break; }
                case 2: {
                    mean=(a[0]+a[1])/2;
                    period=(t[1]-t[0])*2;
                    atime=t[0];
                    a[0]=a[1];
                    a[1]=a[2];
                    t[0]=t[1];
                    t[1]=t[2];
                    j=1;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=0.50;
                        *po++=atime;
                        *po++=period;
                    }
                    break;
                }
                default: {
                    mean=(a[j-1]+a[j-2])/2;
                    period=(t[j-1]-t[j-2])*2;
                    atime=t[j-2];
                    a[j-2]=a[j];
                    t[j-2]=t[j];
                    j=j-2;
                    if (ampl > 0) {
                        *po++=ampl;
                        *po++=mean;
                        *po++=1.00;
                        *po++=atime;
                        *po++=period;
                        cNr1++;
                    }
                    break;
                }
            }
        }
    }
    cNr2 = 1;
    for (index=0; index<j; index++) {
        ampl=fabs(a[index]-a[index+1])/2;
        mean=(a[index]+a[index+1])/2;
        period=(t[index+1]-t[index])*2;
        atime=t[index];
        if (ampl > 0){
            *po++=ampl;
            *po++=mean;
            *po++=0.50;
            *po++=atime;
            *po++=period;
            cNr2++;
        }
    }
//  /* free the memeory !!!*/
//    mxSetN(array_out, tot_num - cNr);
    nout[0] = cNr1;
    nout[1] = cNr2;
}
/* ++++++++++ END RF5 */


/* ++++++++++ BEGIN findcross */
/*
 * findcross.c - 
 *
 *  Returns indices to level v crossings of argument vector
 *
 * 1998 by Per Andreas Brodtkorb. last modified 23.06-98
 */ 
void findcross(double *y, double v, int *ind, int n, int *info) {
    int i,start, ix=0,dcross=0;
    start=0;
    if  ( y[0]< v){
        dcross=-1; /* first is a up-crossing*/ 
    }
    else if ( y[0]> v){
        dcross=1;  /* first is a down-crossing*/ 
    }
    else if  ( y[0]== v){
        /* Find out what type of crossing we have next time.. */
        for (i=1; i<n; i++) {
            start=i;
            if  ( y[i]< v){
                ind[ix] = i-1; /* first crossing is a down crossing*/
                ix++; 
                dcross=-1; /* The next crossing is a up-crossing*/ 
                goto L120;
            }
            else if  ( y[i]> v){
                ind[ix] = i-1; /* first crossing is a up-crossing*/
                ix++; 
                dcross=1;  /*The next crossing is a down-crossing*/ 
                goto L120;
            }
        }
    }
    L120:
    for (i=start; i<n-1; i++) {
        if (( (dcross==-1) && (y[i]<=v) && (y[i+1] > v)  )  || ((dcross==1 ) && (y[i]>=v) && (y[i+1] < v) ) )  { 
      
            ind[ix] = i;
            ix++;
            dcross=-dcross;
        }  
    }
    info[0] = ix;
    return;
}
/* ++++++++++ END findcross */