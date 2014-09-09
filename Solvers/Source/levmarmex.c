/* LEVMARMEX - A MATLAB MEX Interface to LEVMAR
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.i2c2.aut.ac.nz
 */

/* Based in parts on levmar.c supplied with LEVMAR */

#include "mex.h"
#include "mkl.h"
#include <levmar.h>
#include <float.h>

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN], g[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     mxArray *prhs_g[MAXRHS];
     int xrhs, nrhs, xrhs_g, nrhs_g, print;
     double *ydata;
} user_function_data;

//Iteration callback structure
typedef struct {
    char f[FLEN];
    mxArray *plhs[1];
    mxArray *prhs[3];
    bool enabled;
} iter_fun_data;

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }

//Defines
#define MIN_UNCONSTRAINED     0
#define MIN_CONSTRAINED_BC    1
#define MIN_CONSTRAINED_LIC   2
#define MIN_CONSTRAINED_BLIC  3
#define MIN_CONSTRAINED_LEC   4
#define MIN_CONSTRAINED_BLEC  5
#define MIN_CONSTRAINED_LEIC  6
#define MIN_CONSTRAINED_BLEIC 7

//These seem to improve performance on some problems... to be investigated!
#define J_INIT_MU 1e-2
#define J_STOP_THRESH 1e-7

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs, int *conMode);
double getStatus(double stat);
static void func(double *p, double *hx, int m, int n, void *adata);
static void jac(double *p, double *j, int m, int n, void *adata);

//Globals
int citer;
iter_fun_data iterF;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    user_function_data fun;
    double *x0, *ydata = NULL, *lb = NULL, *ub = NULL, *A = NULL, *b = NULL, *Aeq = NULL, *beq = NULL;
    //Options
    int maxIter = 500;
    double info[LM_INFO_SZ];
    double opts[LM_OPTS_SZ]={J_INIT_MU, J_STOP_THRESH, J_STOP_THRESH, J_STOP_THRESH, LM_DIFF_DELTA};
    
    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval;
    double *pcovar = NULL;
    
    //Internal Vars
    size_t ndec, ndat;   
    int i, status, havJac = 0, conMode = 0;
    int nineq=0, neq=0;
    double *covar = NULL;
    double *Apr, *bpr;
    double *llb, *lub;
    citer = 1;
    iterF.enabled = false;
    
    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(LM_VERSION);   
        return;
    }

    //Check user inputs & get constraint information
    checkInputs(prhs,nrhs,&conMode);

    //Get Sizes
    ndec = mxGetNumberOfElements(prhs[2]);
    ndat = mxGetNumberOfElements(prhs[3]);
    //Get Objective Function Handle
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
    fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    fun.print = 0;
    //Check and Get Gradient Function Handle
    if(!mxIsEmpty(prhs[1])) {  
        havJac = 1;
        if (mxIsChar(prhs[1])) {
            CHECK(mxGetString(prhs[1], fun.g, FLEN) == 0,"error reading gradient name string");
            fun.nrhs_g = 1;
            fun.xrhs_g = 0;
        } else {
            fun.prhs_g[0] = (mxArray*)prhs[1];
            strcpy(fun.g, "feval");
            fun.nrhs_g = 2;
            fun.xrhs_g = 1;
        }   
        fun.prhs_g[fun.xrhs_g] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    }

    //Get x0 + data
    x0 = mxGetPr(prhs[2]);
    ydata = mxGetPr(prhs[3]);
    fun.ydata = ydata;
    
    //Get Bounds
    if(conMode & 1) {
        //LB
        if(!mxIsEmpty(prhs[4])){
            llb = mxGetPr(prhs[4]);
            lb = mxCalloc(ndec,sizeof(double));
            memcpy(lb,llb,ndec*sizeof(double));
            for(i=0;i<ndec;i++) {
                if(mxIsInf(lb[i]))
                    lb[i] = -DBL_MAX;
            }
        }
        else {
            lb = mxCalloc(ndec,sizeof(double));
            for(i=0;i<ndec;i++)
                lb[i] = -DBL_MAX;
        }
        //UB
        if(nrhs > 5 && !mxIsEmpty(prhs[5])){
            lub = mxGetPr(prhs[5]);
            ub = mxCalloc(ndec,sizeof(double));
            memcpy(ub,lub,ndec*sizeof(double));
            for(i=0;i<ndec;i++) {
                if(mxIsInf(ub[i]))
                    ub[i] = DBL_MAX;
            }
        }
        else {
            ub = mxCalloc(ndec,sizeof(double));
            for(i=0;i<ndec;i++)
                ub[i] = DBL_MAX;
        }
    }
    //Get Linear Inequality Constraints
    if(conMode & 2) {
        nineq = (int)mxGetM(prhs[7]);
        Apr = mxGetPr(prhs[6]);
        bpr = mxGetPr(prhs[7]);
        //Need to flip >= to <=
        A = mxCalloc(ndec*nineq,sizeof(double));
        b = mxCalloc(nineq,sizeof(double));
        for(i=0;i<ndec*nineq;i++)
            A[i] = -Apr[i];
        for(i=0;i<nineq;i++)
            b[i] = -bpr[i];
    }
    //Get Linear Equality Constraints
    if(conMode & 4) {
        Aeq = mxGetPr(prhs[8]);
        beq = mxGetPr(prhs[9]);
        neq = (int)mxGetM(prhs[9]);
    }
    
    //Get Options if specified
    if(nrhs > 10) {
        if(mxGetField(prhs[10],0,"maxiter"))
            maxIter = (int)*mxGetPr(mxGetField(prhs[10],0,"maxiter"));
        if(mxGetField(prhs[10],0,"display"))
            fun.print = (int)*mxGetPr(mxGetField(prhs[10],0,"display"));
        if(mxGetField(prhs[10],0,"iterfun") && !mxIsEmpty(mxGetField(prhs[10],0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(prhs[10],0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
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
    
    //Copy initial guess to x
    memcpy(x,x0,ndec*sizeof(double));
    //Create Covariance Matrix if Required
    if(nlhs>4)
        covar=mxCalloc(ndec*ndec,sizeof(double));
    
    //Print Header
    if(fun.print) {
        mexPrintf("\n------------------------------------------------------------------\n");        
        mexPrintf(" This is LEVMAR v%s\n",LM_VERSION);           
        mexPrintf(" Author: Manolis Lourakis\n\n");
        mexPrintf(" Algorithm: ");
        switch(conMode)
        {
            case MIN_UNCONSTRAINED: mexPrintf("Unconstrained"); break;
            case MIN_CONSTRAINED_BC: mexPrintf("Box Constrained"); break;
            case MIN_CONSTRAINED_LIC: mexPrintf("Linear Inequality"); break;
            case MIN_CONSTRAINED_BLIC: mexPrintf("Boxed Linear Inequality"); break;
            case MIN_CONSTRAINED_LEC: mexPrintf("Linear Equality"); break;
            case MIN_CONSTRAINED_BLEC: mexPrintf("Boxed Linear Equality"); break;
            case MIN_CONSTRAINED_LEIC: mexPrintf("Linear Inequality + Equality"); break;
            case MIN_CONSTRAINED_BLEIC: mexPrintf("Boxed Linear Inequality + Equality"); break;
        }
        if(havJac)
            mexPrintf(" (with Jacobian)\n\n");
        else
            mexPrintf("\n\n");
        
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Data Points:            %4d\n",ndat);
        mexPrintf(" # Inequality Constraints: %4d\n",nineq);
        mexPrintf(" # Equality Constraints:   %4d\n",neq);
        mexPrintf("------------------------------------------------------------------\n");
    }
  
    //Solve based on constraints
    switch(conMode)
    {
        case MIN_UNCONSTRAINED:
            if(havJac)
                status = dlevmar_der(func, jac, x, ydata, (int)ndec, (int)ndat, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_dif(func, x, ydata, (int)ndec, (int)ndat, maxIter, opts, info, NULL, covar, &fun);            
            break;
        case MIN_CONSTRAINED_BC:
            if(havJac)
                status = dlevmar_bc_der(func, jac, x, ydata, (int)ndec, (int)ndat, lb, ub, NULL, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_bc_dif(func, x, ydata, (int)ndec, (int)ndat, lb, ub, NULL, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_LIC:
            if(havJac)
                status = dlevmar_lic_der(func, jac, x, ydata, (int)ndec, (int)ndat, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_lic_dif(func, x, ydata, (int)ndec, (int)ndat, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_BLIC:
            if(havJac)
                status = dlevmar_blic_der(func, jac, x, ydata, (int)ndec, (int)ndat, lb, ub, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_blic_dif(func, x, ydata, (int)ndec, (int)ndat, lb, ub, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_LEC:
            if(havJac)
                status = dlevmar_lec_der(func, jac, x, ydata, (int)ndec, (int)ndat, Aeq, beq, neq, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_lec_dif(func, x, ydata, (int)ndec, (int)ndat, Aeq, beq, neq, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_BLEC:
            if(havJac)
                status = dlevmar_blec_der(func, jac, x, ydata, (int)ndec, (int)ndat, lb, ub, Aeq, beq, neq, NULL, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_blec_dif(func, x, ydata, (int)ndec, (int)ndat, lb, ub, Aeq, beq, neq, NULL, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_LEIC:
            if(havJac)
                status = dlevmar_leic_der(func, jac, x, ydata, (int)ndec, (int)ndat, Aeq, beq, neq, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_leic_dif(func, x, ydata, (int)ndec, (int)ndat, Aeq, beq, neq, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            break;
        case MIN_CONSTRAINED_BLEIC:
            if(havJac)
                status = dlevmar_bleic_der(func, jac, x, ydata, (int)ndec, (int)ndat, lb, ub, Aeq, beq, neq, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            else
                status = dlevmar_bleic_dif(func, x, ydata, (int)ndec, (int)ndat, lb, ub, Aeq, beq, neq, A, b, nineq, maxIter, opts, info, NULL, covar, &fun);
            break;
        default:
            mexErrMsgTxt("Unknown constraint configuration");
    }
       
    //Save Status & Iterations
    *fval = info[1];
    *exitflag = info[6];
    *iter = (double)status;
    *feval = (double)citer;
    
    //Save Covariance if Required
    if(nlhs > 5) {
        plhs[5] = mxCreateDoubleMatrix(ndec, ndec, mxREAL);
        pcovar = mxGetPr(plhs[5]);
        memcpy(pcovar,covar,ndec*ndec*sizeof(double));
    }
    
    //Print Header
    if(fun.print){ 
        switch((int)*exitflag)
        {
            case 1: mexPrintf("\n *** SUCCESSFUL TERMINATION: Gradient Convergence ***\n"); break;
            case 2: mexPrintf("\n *** SUCCESSFUL TERMINATION: Dp Convergence ***\n"); break;
            case 3: mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;
            case 4: mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Singular Matrix (Try Increasing Mu) ***\n"); break;
            case 5: mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: No Further Error Reduction Possible (Try Increasing Mu) ***\n"); break;
            case 6: mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Stopped By Small ||e||_2 ***\n"); break;
            case 7: mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Detected NaN or Inf values ***\n"); break;
            default: mexPrintf("\n *** TERMINATION: ROUTINE ERROR (Code %d) ***\n",(int)*exitflag); break;
        }
        if(*exitflag==1)
            mexPrintf(" Final SSE: %12.5g\n In %3.0f iterations\n",*fval,*iter);
        
        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Clean Up
    if(lb) mxFree(lb);
    if(ub) mxFree(ub);
    if(covar) mxFree(covar);
    if(A) mxFree(A);
    if(b) mxFree(b);
}

static void func(double *p, double *hx, int m, int n, void *adata)
{
    bool havrnorm = false, stop = false;
    int stat, i;
    double *fval, rnorm;
    user_function_data *fun = (user_function_data *) adata;

    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), p, m * sizeof(double));

    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");

    //Get Objective
    fval = mxGetPr(fun->plhs[0]);
    memcpy(hx,fval,n*sizeof(double));
    
    //Iteration Printing
    if(fun->print > 1) {
        if(citer == 1 || !(citer%10))
            mexPrintf(" feval             sse\n");
        
        //Calculate Residual Norm
        rnorm = 0; havrnorm = true;
        for(i=0;i<n;i++)
            rnorm += (fval[i]-fun->ydata[i])*(fval[i]-fun->ydata[i]);

        mexPrintf("%5d     %12.5g\n",citer,rnorm);        
        mexEvalString("drawnow;"); //flush draw buffer  
    }
    
    //Iteration Callback
    if(iterF.enabled)
    {
        //Calculate sse if we don't have it
        if(!havrnorm)
            for(i=0;i<n;i++)
                rnorm += (fval[i]-fun->ydata[i])*(fval[i]-fun->ydata[i]);

        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), &citer, sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), &rnorm, sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), p, m * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
    }
    
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
        
    citer++;
    
    //Warn user stop not implemented
    if(stop)
        mexWarnMsgTxt("LEVMAR does not implement the stop feature of iterfun");
}

static void jac(double *p, double *j, int m, int n, void *adata)
{
    int i, k, stat;
    double *grad;
    user_function_data *fun = (user_function_data *) adata;

    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs_g[fun->xrhs_g]), p, m * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs_g, fun->prhs_g, fun->g);
    if(stat)
      mexErrMsgTxt("Error calling Gradient Function!");

    //Get Gradient & Transpose
    grad = mxGetPr(fun->plhs[0]);
    for(i=0; i<n; ++i)
        for(k=0; k<m; ++k)
            j[i*m+k]=grad[i+k*n];

    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
}

void checkInputs(const mxArray *prhs[], int nrhs, int *conMode)
{    
    size_t Mx0;
    
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 4 arguments to levmar!\n\nlevmar(fun,grad,x0,ydata) or\nlevmar(fun,grad,x0,ydata,lb,ub,A,b,Aeq,beq,opts)");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(prhs[0]))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsEmpty(prhs[1]) && (!mxIsFunctionHandle(prhs[1]) && !mxIsChar(prhs[1])))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsEmpty(prhs[2]))
        mexErrMsgTxt("x0 must be a real double column vector!");
    if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsEmpty(prhs[3]))
        mexErrMsgTxt("ydata must be a real double column vector!");
    
    //Get ndec
    Mx0 = mxGetNumberOfElements(prhs[2]);
    
    //Check Bounds
    if(nrhs > 4) {
        if(!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
            mexErrMsgTxt("lb must be a real double column vector!");
        if(nrhs > 5 && (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])))
            mexErrMsgTxt("ub must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(prhs[4]) && (Mx0 != mxGetNumberOfElements(prhs[4])))
            mexErrMsgTxt("lb is not the same length as x0! Ensure they are both Column Vectors");
        if(nrhs > 5 && !mxIsEmpty(prhs[5]) && (Mx0 != mxGetNumberOfElements(prhs[5])))
            mexErrMsgTxt("ub is not the same length as x0! Ensure they are both Column Vectors");
        //Confirm Bounds
        if(!mxIsEmpty(prhs[4]) || (nrhs > 5 && !mxIsEmpty(prhs[5])))
            *conMode |= 1;
    }
    
    //Check Linear Inequality
    if(nrhs > 6) {
        if(nrhs == 7)
            mexErrMsgTxt("You must supply both A & b!");
        if(!mxIsEmpty(prhs[6]) && (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || mxIsSparse(prhs[6])))
            mexErrMsgTxt("A must be a full, real double matrix!");
        if(!mxIsEmpty(prhs[7]) && (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7])))
            mexErrMsgTxt("b must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(prhs[6])) {
            if(mxGetM(prhs[6]) != Mx0)
                mexErrMsgTxt("A is not the right size! Remember to transpose it for levmar.");
            if(mxGetN(prhs[6]) != mxGetM(prhs[7]))
                mexErrMsgTxt("A & b sizes do not correspond");
            //Confirm Lin Eq
            *conMode |= 2;  
        }
    }    
    
    //Check Linear Equality
    if(nrhs > 8) {
        if(nrhs == 9)
            mexErrMsgTxt("You must supply both Aeq & beq!");
        if(!mxIsEmpty(prhs[8]) && (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8]) || mxIsSparse(prhs[8])))
            mexErrMsgTxt("Aeq must be a full, real double matrix!");
        if(!mxIsEmpty(prhs[9]) && (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9])))
            mexErrMsgTxt("beq must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(prhs[8])) {
            if(mxGetM(prhs[8]) != Mx0)
                mexErrMsgTxt("Aeq is not the right size! Remember to transpose it for levmar.");
            if(mxGetN(prhs[8]) != mxGetM(prhs[9]))
                mexErrMsgTxt("Aeq & beq sizes do not correspond");
            //Confirm Lin Eq
            *conMode |= 4;    
        }
    }
    
    //Check Options
    if(nrhs > 10) {
        if(!mxIsStruct(prhs[10]))
            mexErrMsgTxt("The specified options must be a structure!");
    }

}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" LEVMAR: Levenberg-Marquardt Nonlinear Least Squares in C/C++ [v%s, Built %s]\n",LM_VERSION,__DATE__);
    mexPrintf("  - Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html\n");
    mexPrintf("  - Source available from: http://www.ics.forth.gr/~lourakis/levmar/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    
    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}