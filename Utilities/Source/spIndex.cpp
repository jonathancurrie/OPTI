#include "mex.h"
#include <stdio.h>
#include <string.h>
#define max(a, b)  (((a) > (b)) ? (a) : (b))

void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] ) 
{
    if(nrhs < 1 || !mxIsSparse(prhs[0]))
        mexErrMsgTxt("This function requires one sparse input");

    mwIndex *ir = mxGetIr(prhs[0]);
    mwIndex *jc = mxGetJc(prhs[0]);
    double *pr = mxGetPr(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    mwSize m = mxGetM(prhs[0]);
    mwSize nnz = jc[n];
    
    if(nlhs >= 3) {
        plhs[0] = mxCreateDoubleMatrix(n+1,1,mxREAL);
        plhs[1] = mxCreateDoubleMatrix(nnz,1,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(nnz,1,mxREAL);
        double *jc_out = mxGetPr(plhs[0]);
        double *ir_out = mxGetPr(plhs[1]);
        double *pr_out = mxGetPr(plhs[2]);
        
        for(int i = 0; i <= n; i++)
            jc_out[i] = (double)jc[i];
        for(int i = 0; i < nnz; i++)
            ir_out[i] = (double)ir[i];                  
        memcpy(pr_out,pr,nnz*sizeof(double));       
    }
    else {
        if(n+1 > nnz) {
            for(int i = 0; i < n+1; i++) {
                if(i < nnz)
                    mexPrintf("jc[%d] %3d, ir[%d] %3d, pr[%d] %f\n",i,jc[i],i,ir[i],i,pr[i]);
                else
                    mexPrintf("jc[%d] %3d\n",i,jc[i]);

            }
        } else {
            for(int i = 0; i < nnz; i++) {
                if(i <= n)
                    mexPrintf("jc[%d] %3d, ir[%d] %3d, pr[%d] %f\n",i,jc[i],i,ir[i],i,pr[i]);
                else
                    mexPrintf("           ir[%d] %3d, pr[%d] %f\n",i,ir[i],i,pr[i]);

            }
        }            
    }
}