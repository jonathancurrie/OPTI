#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <csignal>
#include "lbfgsb_program.h"

// Function definitions.
// -----------------------------------------------------------------
int getBoundType (double& lb, double& ub) {
  bool haslb = !mxIsInf(lb);
  bool hasub = !mxIsInf(ub);
  int  btype = haslb + 3*hasub - 2*(haslb && hasub);

  lb = haslb ? lb : 0;
  ub = hasub ? ub : 0;

  return btype;
}

// Function definitions for class L-BFGS-B Program.
// -----------------------------------------------------------------
lbfgsb_program::lbfgsb_program (user_function_data *fun, user_function_data *grad, iter_fun_data *iterF, size_t ndec, 
                                double *in_lb, double *in_ub, double *in_x0, double *xfinal, 
                                double *fval, double *iter, int printLevel, int maxIter, int maxFeval, double maxtime, double ftol, int m, double pgtol)

    : Program((int)ndec,0,0,0,0,m,maxIter,maxFeval,maxtime,ftol/mxGetEps(),pgtol) //Program Constructor
{
    //Create local copies of input args (we will change them)
    x     = new double[ndec];
    lb    = new double[ndec];
    ub    = new double[ndec];
    for(int i = 0; i < ndec; i++) {
        lb[i] = in_lb[i];
        ub[i] = in_ub[i];
        x[i] = in_x0[i];
    }
    //Solver bound types
    btype = new int[ndec];
    // Set the bound types.
    for (int i = 0; i < ndec; i++)
        btype[i] = getBoundType(lb[i],ub[i]);
    //Save size
    this->ndec = ndec;
    //Save print stuff
    this->printLevel = printLevel;
    
    //Save Function Information
    this->fun = fun;
    this->grad = grad;
    this->iterF = iterF;
    
    //Assign Output Variables
    this->xfinal = xfinal;

    this->noFevals = 0;
}

lbfgsb_program::~lbfgsb_program() {
    delete[] x;
    delete[] lb;
    delete[] ub;
    delete[] btype;
}

double lbfgsb_program::computeObjective (int n, double* x) 
{     
    double f = 0;

    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, n * sizeof(double));
    noFevals++;
    
    try {
    	mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
        //Get Objective
        f = mxGetScalar(fun->plhs[0]);
        // Clean up Ptr
        mxDestroyArray(fun->plhs[0]);
    }
    catch(...)
    {
        mexWarnMsgTxt("Unrecoverable Error from Objective Callback, Exiting LBFGSB...\n");
        //Force exit
        utSetInterruptPending(true);
    }
    return f;
}

void lbfgsb_program::computeGradient (int n, double* x, double* g) 
{
    double *gval;

    grad->plhs[0] = NULL;
    memcpy(mxGetPr(grad->prhs[grad->xrhs]), x, n * sizeof(double));

    try {
    	mexCallMATLAB(1, grad->plhs, grad->nrhs, grad->prhs, grad->f);
        //Get Gradient
        gval = mxGetPr(grad->plhs[0]);
        memcpy(g,gval,n*sizeof(double));

        // Clean up Ptr
        mxDestroyArray(grad->plhs[0]);
    }
    catch(...)
    {
        mexWarnMsgTxt("Unrecoverable Error from Gradient Callback, Exiting LBFGSB...\n");
        //Force exit
        utSetInterruptPending(true);
    } 
}

bool lbfgsb_program::iterCallback (int t, double etime, double* x, double f) 
{
    bool stop = false;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting L-BFGS-B...\n\n");
        return true; //terminate asap
    }
    
    //Iteration Printing
    if(printLevel > 1) {      
        if(t==1 || !(t%10))
            mexPrintf(" iter      feval     time[s]      objective\n");
        
        mexPrintf("%5d     %5d       %5.2f    %12.5g\n",t,noFevals,etime,f);       
        mexEvalString("drawnow;"); //flush draw buffer        
    }
    
    //Iteration Callback
    if(iterF->enabled)
    {
        iterF->plhs[0] = NULL;
        memcpy(mxGetData(iterF->prhs[1]), &t, sizeof(int));
        memcpy(mxGetPr(iterF->prhs[2]), &f, sizeof(double));
        memcpy(mxGetPr(iterF->prhs[3]), x, this->ndec * sizeof(double));
        try {
        	mexCallMATLAB(1, iterF->plhs, 4, iterF->prhs, iterF->f);
        }
        catch(...)
        {
            mexWarnMsgTxt("Unrecoverable Error from Iteration Callback, Exiting LBFGSB...\n");
            //Force exit
            return true;
        }

        //Collect return argument
        stop = *(bool*)mxGetData(iterF->plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF->plhs[0]);
    }
    
    return stop;
}

int lbfgsb_program::runSolver(int &iter, int &feval, double &fval) 
{
    SolverExitStatus status;  // The solver success.
    double evaltime;

    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is L-BFGS-B v3.0 (August 2011)\n Authors: C. Zhu, R. Byrd, J. Nocedal, J. L. Morales\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);

        if(printLevel > 1)
            mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;");
    }
  
    status = Program::runSolver(iter,fval,evaltime);

    // Copy the solution
    for(int i = 0; i < ndec; i++)
      xfinal[i] = x[i];

    //Get no Fevals
    feval = noFevals;
    
    //Print Footer
    if(printLevel > 0){            
        //Termination Detected
        switch(status)
        {
            case success: mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n"); break;
            case exceededIterations: mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;
            case exceededMaxFeval: mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); break;
            case exceededMaxTime: mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;
            case abnormalTermination: mexPrintf("\n *** TERMINATION: ABNORMAL TERMINATION ***\n"); break;
            case errorOnInput: mexPrintf("\n *** TERMINATION: ERROR ON INPUT ***\n"); break;
            case userExit: mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
            default: mexPrintf("\n *** TERMINATION: STATUS UNKNOWN ***\n"); break;
        }        
        mexPrintf(" Objective Value: %1.6g\n",fval);
        mexPrintf(" In: %d iterations and %1.2f seconds\n",iter,evaltime);
        mexPrintf("------------------------------------------------------------------\n\n");
    }

    return status;
}
