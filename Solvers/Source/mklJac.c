/* MKLJAC - A simple MEX Interface to the MKL djacobi Function
 * Copyright (C) 2011 Jonathan Currie (IPL)                             
 */

#include <mex.h>
#include <mkl.h>
#include <string.h>

//Matlab Data Structure
typedef struct {
     mxArray *plhs[1];
     mxArray *prhs[2];
} matlab_data;

void checkInputs(int nrhs, const mxArray *prhs[]);
MKL_INT dummyCallForSize(const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *x;
    MKL_INT n, m = 1;
    
    //Return Args
    double *fjac, *status;
    
    //Internal Args
    extern void mFun(MKL_INT*, MKL_INT*, double*, double*, void*);
    double tol = 1e-6;
    matlab_data mdata;
    
    //Check Inputs
    if(nrhs < 1) {
        mexPrintf("This is a MEX interface to Intel MKL %d.%d R%d djacobi function\n\nUsage: [jac,status] = mklJac(fun,x,nrow)\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
        return;
    }
    checkInputs(nrhs,prhs);    
    
    //Get Inputs
    x = mxGetPr(prhs[1]); 
    n = (MKL_INT)mxGetNumberOfElements(prhs[1]);       //length of x
    if(nrhs > 2 && !mxIsEmpty(prhs[2]))
        m = (MKL_INT)*mxGetPr(prhs[2]); //length of f 
    else
        m = dummyCallForSize(prhs);     //have to get length from a dummy call (bad...)
    if(nrhs > 3 && !mxIsEmpty(prhs[3]))
        tol = *mxGetPr(prhs[3]);
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL); 
    fjac = mxGetPr(plhs[0]); 
    status = mxGetPr(plhs[1]); 

    //Build Data Structure
    mdata.prhs[0] = prhs[0];
    mdata.prhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);

    //Call Jacobi Function
    if(djacobix(mFun,&n,&m,fjac,x,&tol,&mdata) == TR_SUCCESS)
        *status = 1;
    else
        *status = 0;
}

//Wrapper to Call Matlab Function
void mFun(MKL_INT *m, MKL_INT *n, double *x, double *f, void *data)
{
    matlab_data *mdata = (matlab_data*)data;

    //Clear Current Pointer
    mdata->plhs[0] = NULL;
    //Insert x into RHS
    memcpy(mxGetPr(mdata->prhs[1]), x, *n * sizeof(double));
    //Call Matlab Function
    mexCallMATLAB(1, mdata->plhs, 2, mdata->prhs, "feval");
    //Check Return Size & Data Type
    if(mxGetNumberOfElements(mdata->plhs[0]) != *m)
        mexErrMsgTxt("Incorrect sized vector returned from function - check nrow is specified correctly");
    if(mxIsSparse(mdata->plhs[0]) || !mxIsDouble(mdata->plhs[0]) || mxIsComplex(mdata->plhs[0]))
        mexErrMsgTxt("Your function must only return real, dense, double precision numeric data!");
    //Save Result
    memcpy(f, mxGetPr(mdata->plhs[0]), *m * sizeof(double));    
    //Clean up
    mxDestroyArray(mdata->plhs[0]);
}

void checkInputs(int nrhs, const mxArray *prhs[])
{
    MKL_INT m;
    double tol;
    
    if(nrhs < 2)
        mexErrMsgTxt("You must supply at least 2 arguments to mklJac!");
    
    //fun
    if(mxIsEmpty(prhs[0]))
        mexErrMsgTxt("fun is empty!");
    if(!mxIsFunctionHandle(prhs[0]))
        mexErrMsgTxt("fun must be a function handle!");
    //x
    if(mxIsEmpty(prhs[1]))
        mexErrMsgTxt("x is empty!");
    //m
    if(nrhs > 2 && !mxIsEmpty(prhs[2])) {
        m = (MKL_INT)*mxGetPr(prhs[2]);
        if(m < 1 || m > 1e8)
            mexErrMsgTxt("The value of nrow must be 0 < len < 1e8");
    }
    //tol
    if(nrhs > 3 && !mxIsEmpty(prhs[3])) {
        tol = *mxGetPr(prhs[3]);
        if(tol < 1e-15 || tol > 1)
            mexErrMsgTxt("tol must be 1e-15 < tol < 1");
    }
}

//Dummy call to get size
MKL_INT dummyCallForSize(const mxArray *prhs[])
{
    MKL_INT m;
    matlab_data mdata;
    
    mdata.plhs[0] = NULL;
    mdata.prhs[0] = prhs[0];
    mdata.prhs[1] = prhs[1];
    mexCallMATLAB(1, mdata.plhs, 2, mdata.prhs, "feval");
    m = (MKL_INT)mxGetNumberOfElements(mdata.plhs[0]);
    mxDestroyArray(mdata.plhs[0]);
    return m;
}
