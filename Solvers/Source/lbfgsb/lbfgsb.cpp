/* LBFGSB - A simple MEX Interface to L-BFGS-B NLP Solver
 * Copyright (C) Jonathan Currie 2011 (IPL)                             
 */

// Based largely on parts by Dr. Peter Carbonetto

#include "mex.h"
#include "lbfgsb.h"
#include "lbfgsb_program.h"
#include <exception>
#include "opti_util.h"
#ifdef LINK_MKL
    #include "mkl.h"
#endif
#include <string>

#define LBFGSB_VERSION "3.0"

//Function Prototypes
void printSolverInfo();
void GetIntegerOption(const mxArray *opts, char *name, int *var);
void GetDoubleOption(const mxArray *opts, char *name, double *var);
void checkInputs(const mxArray *prhs[], int nrhs);

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    user_function_data fun, grad;
    iter_fun_data iterF;
    double *x0, *lb, *ub;
    //Options
    int printLevel = 0, maxIter = 1000, nupdate = 5, maxFeval = 1500;    
    double ftol = 1e-7, pgtol = 1e-5, maxtime = 1000;
    iterF.enabled = false;
    
    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval;
    
    //Internal Vars
    size_t ndec;
    int iiter, nfeval;
    
    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString(LBFGSB_VERSION);
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs);

    try 
    {
        //Get Size
        ndec = mxGetNumberOfElements(prhs[4]);
        //Get Function Handles
        if (mxIsChar(prhs[0])) {
            CHECK(mxGetString(prhs[0], fun.f, FLEN) == 0,"error reading objective name string");
            fun.nrhs = 1;
            fun.xrhs = 0;
        } else {
            fun.prhs[0] = (mxArray*)prhs[0];
            strcpy(fun.f, "feval");
            fun.nrhs = 2;
            fun.xrhs = 1;
        }
        if (mxIsChar(prhs[1])) {
            CHECK(mxGetString(prhs[1], grad.f, FLEN) == 0,"error reading gradient name string");
            grad.nrhs = 1;
            grad.xrhs = 0;
        } else {
            grad.prhs[0] = (mxArray*)prhs[1];
            strcpy(grad.f, "feval");
            grad.nrhs = 2;
            grad.xrhs = 1;
        }
        fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
        grad.prhs[grad.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL);
        
        //Get Bounds + x0
        lb = mxGetPr(prhs[2]);
        ub = mxGetPr(prhs[3]);
        x0 = mxGetPr(prhs[4]);
        
        //Get Options if specified
        if(nrhs > 5) {
            GetIntegerOption(prhs[5], "display", &printLevel);
            GetIntegerOption(prhs[5], "maxiter", &maxIter);
            GetIntegerOption(prhs[5], "maxfeval", &maxFeval);
            GetIntegerOption(prhs[5], "nupdate", &nupdate);
            GetDoubleOption(prhs[5], "tolrfun", &ftol);
            GetDoubleOption(prhs[5], "tolrfun", &ftol);
            GetDoubleOption(prhs[5], "pgtol", &pgtol);
            GetDoubleOption(prhs[5], "maxtime", &maxtime);            
            if(mxGetField(prhs[5],0,"iterfun") && !mxIsEmpty(mxGetField(prhs[5],0,"iterfun")))
            {
                iterF.prhs[0] = (mxArray*)mxGetField(prhs[5],0,"iterfun");
                strcpy(iterF.f, "feval");
                iterF.enabled = true;  
                iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
                iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
                iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
            }
            CheckOptiVersion(prhs[5]);
        }                        
        
        //Create Outputs
        plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
        plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
        x = mxGetPr(plhs[0]); 
        fval = mxGetPr(plhs[1]); 
        exitflag = mxGetPr(plhs[2]);    
        iter = mxGetPr(plhs[3]);
        feval = mxGetPr(plhs[4]);
        
        //Create Class for Peter's L-BFGS-B Interface
        lbfgsb_program solver(&fun,&grad,&iterF,ndec,lb,ub,x0,x,fval,iter,printLevel,maxIter,maxFeval,maxtime,ftol,nupdate,pgtol);
        
        //Run the Solver
        int exitStatus = solver.runSolver(iiter,nfeval,*fval);
        //Save Status & Iterations
        *exitflag = (double)exitStatus;
        *iter = (double)iiter;
        *feval = (double)nfeval;

    } 
    catch (std::exception& error) 
    {
        mexErrMsgTxt(error.what());
    }
}

void checkInputs(const mxArray *prhs[], int nrhs)
{    
    size_t Mlb, Mub, Mx0;
    
    if(nrhs < 5)
        mexErrMsgTxt("You must supply at least 5 arguments to lbfsgb!\nlbfgsb(fun,grad,lb,ub,x0)");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(prhs[0]))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsFunctionHandle(prhs[1]) && !mxIsChar(prhs[1]))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsEmpty(prhs[2]))
        mexErrMsgTxt("lb must be a real double column vector!");
    if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsEmpty(prhs[3]))
        mexErrMsgTxt("ub must be a real double column vector!");
    if(!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxIsEmpty(prhs[4]))
        mexErrMsgTxt("x0 must be a real double column vector!");
    if((nrhs > 5) && !mxIsStruct(prhs[5]))
        mexErrMsgTxt("Options must be supplied as a structure!");
    
    //Check Sizes
    Mlb = mxGetNumberOfElements(prhs[2]);
    Mub = mxGetNumberOfElements(prhs[3]);
    Mx0 = mxGetNumberOfElements(prhs[4]);
    
    if(Mlb != Mub)
        mexErrMsgTxt("The bound vectors are not the same length!");
    if(Mlb != Mx0)
        mexErrMsgTxt("x0 is not the same length as the bounds!");
}

//Print Solver Information
void printSolverInfo()
{     
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" L-BFGS-B: Limited Memory Broyden-Fletcher-Goldfarb-Shanno Bounded Optimization [%s]\n",LBFGSB_VERSION);
    PRINT_BUILD_INFO;
    mexPrintf("  - Released under the BSD 3 Clause License: http://en.wikipedia.org/wiki/BSD_licenses\n");
    mexPrintf("  - Source available from: http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    #endif
    #ifdef LINK_NETLIB_BLAS
        mexPrintf("  - NETLIB BLAS: http://www.netlib.org/blas/\n  - NETLIB LAPACK: http://www.netlib.org/lapack/\n");
    #endif
    
    mexPrintf("\n MEX Interface P.Carbonetto [Modified by J.Currie 2013]\n");
    mexPrintf("-----------------------------------------------------------\n");
}

//Option Getting Methods
void GetIntegerOption(const mxArray *opts, char *name, int *var)
{
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
        *var = (int)*mxGetPr(mxGetField(opts,0,name));
}
void GetDoubleOption(const mxArray *opts, char *name, double *var)
{
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
        *var = *mxGetPr(mxGetField(opts,0,name));
}