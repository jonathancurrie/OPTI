/* LMDEREX - A MATLAB MEX Interface to LMDER
 * Released Under the BSD 3-Clause License:
 * https://www.controlengineering.co.nz/Wikis/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2011-2013
 * www.inverseproblem.co.nz
 */

#include "mex.h"
#include <time.h>
#include "opti_build_utils.h"
#ifdef LINK_MKL
    #include "mkl.h"
#endif

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     int xrhs, nrhs;
} user_function_data;

//Iteration callback structure
typedef struct {
    char f[FLEN];
    mxArray *plhs[1];
    mxArray *prhs[3];
    bool enabled;
} iter_fun_data;

//Ctrl-C Detection
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
double getStatus(int stat);
static void func(int *m, int *n, double *x, double *fvec, int *iflag);
static void jacfunc(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag);

// LMDIF routine implemented in Fortran 77.
extern void LMDIF (void (*fun)(int*,int*,double*,double*,int*),int *m, int *n, double* x, double *fvec, double *ftol, 
                   double *xtol, double *gtol, int *maxfev, double *espfcn, double *diag, int *mode, double *factor, 
                   int *nprint, int *info, int *nfev, double *fjac, int *ldfjac, int *ipvt, double *qtf, double *wa1, 
                   double *wa2, double *wa3, double *wa4);

// LMDER routine implemented in Fortran 77.
extern void LMDER (void (*fun)(int*,int*,double*,double*,double*,int*,int*),int *m, int *n, double* x, double *fvec, 
                   double *fjac, int *ldfjac, double *ftol, double *xtol, double *gtol, int *maxfev, double *diag, 
                   int *mode, double *factor, int *nprint, int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
                   double *wa1, double *wa2, double *wa3, double *wa4);


//Global Variables (no user arg passing to LMDIF functions...)
user_function_data fun;
user_function_data grad;
iter_fun_data iterF;
double *ydata;
int printLevel, citer;
double maxtime;
clock_t start, end;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Outputs Args
    double *x, *fval, *exitflag, *iter;
    
    //Internal Vars
    double *mptr;
    size_t ndec, ndat;   
    int i, havJac = 0;
    
    //LMDIF + LMDER Vars
    int m;                  //len data
    int n;                  //len x
    double *x0;             //initial guess + sol
    double *fvec;           //final fvec
    double ftol = 1e-7;     //function tolerance
    double xtol = 1e-7;     //iterate tolerance
    double gtol = 1e-7;     //gradient tolerance
    int maxfev = 1500;      //max fevals
    double epsfcn = mxGetEps(); //FD max error
    double *diag;           //scaling (ignored in mode 1)   
    int mode = 1;           //internal scaling
    double factor = 100;    //used for initial step bound
    int nprint = -1;        //no printing
    int info = 0;           //output status
    int nfev = 0;           //number of fevals
    int njev = 0;           //number of jevals
    double *fjac;           //final jacobian
    int ldfjac;             //ld of jac
    int *ipvt;              //permutation
    double *qtf;            //qt * fvec    
    double *wa1, *wa2, *wa3, *wa4;  //work vectors

    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString("");  //no version info? 
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs);
    
    //Set Defaults
    citer = 1;
    printLevel = 0;
    iterF.enabled = false;
    maxtime = 1000;

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
    //Check and Get Gradient Function Handle
    if(!mxIsEmpty(prhs[1])) {  
        havJac = 1;
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
        grad.prhs[grad.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    }

    //Get x0 + data
    x0 = mxGetPr(prhs[2]);
    ydata = mxGetPr(prhs[3]);        
    
    //Get Options if specified
    if(nrhs > 4) {
        if(mxGetField(prhs[4],0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(prhs[4],0,"display"));
        if(mxGetField(prhs[4],0,"maxfeval"))
            maxfev = (int)*mxGetPr(mxGetField(prhs[4],0,"maxfeval"));
        if(mxGetField(prhs[4],0,"maxtime"))
            maxtime = *mxGetPr(mxGetField(prhs[4],0,"maxtime"));
        if(mxGetField(prhs[4],0,"tolrfun"))
            ftol = *mxGetPr(mxGetField(prhs[4],0,"tolrfun"));
        if(mxGetField(prhs[4],0,"iterfun") && !mxIsEmpty(mxGetField(prhs[4],0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(prhs[4],0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
        CheckOptiVersion(prhs[4]);
    }       
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    iter = mxGetPr(plhs[3]);
    
    //Copy initial guess to x
    memcpy(x,x0,ndec*sizeof(double));
    
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        if(havJac)
            mexPrintf(" This is LM_DER\n");
        else
            mexPrintf(" This is LM_DIF\n");
            
        mexPrintf(" Authors: Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Data Points:            %4d\n",ndat);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Assign Arguments
    m = (int)ndat; ldfjac = m;
    n = (int)ndec;
    diag = mxCalloc(ndec,sizeof(double));
    fvec = mxCalloc(ndat,sizeof(double));
    fjac = mxCalloc(ndec*ndat,sizeof(double));
    ipvt = mxCalloc(ndec,sizeof(int));
    qtf = mxCalloc(ndec,sizeof(double));   
    //Get Work Memory
    mptr = mxCalloc(3*ndec+ndat,sizeof(double)); i = 0;
    wa1 = &mptr[i]; i += (int)ndec;
    wa2 = &mptr[i]; i += (int)ndec;
    wa3 = &mptr[i]; i += (int)ndec;
    wa4 = &mptr[i];
    
    //Start Timer
    start = clock();
    
    //Call Algorithm based on Derivatives
    if(havJac)
        LMDER(jacfunc,&m,&n,x,fvec,fjac,&ldfjac,&ftol,&xtol,&gtol,&maxfev,diag,&mode,&factor,&nprint,&info,&nfev,&njev,ipvt,qtf,wa1,wa2,wa3,wa4);
    else
        LMDIF(func,&m,&n,x,fvec,&ftol,&xtol,&gtol,&maxfev,&epsfcn,diag,&mode,&factor,&nprint,&info,&nfev,fjac,&ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4);
    
    //Calculate final Rnorm
    *fval = 0;
    for(i=0;i<m;i++)
        *fval += fvec[i]*fvec[i];
    
    //Check if maxtime exceeded
    if(((double)(end-start))/CLOCKS_PER_SEC > maxtime)
        info = 15;
    
    //Save Status & Iterations
    *exitflag = getStatus(info);
    *iter = (double)nfev;
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch((int)info)
        {
            //Success
            case 1:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Both actual and predicted relative reductions in the sum of squares are at most ftol ***\n"); break;
            case 2:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Relative error between two consecutive iterates is at most xtol ***\n"); break;
            case 3:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Reductions in the sum of squares are at most ftol AND relative error between two consecutive iterates is at most xtol ***\n"); break;
            case 4:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value ***\n"); break;
            //Error
            case 0:
                mexPrintf("\n *** ERROR: IMPROPER INPUT PARAMETERS ***\n"); break;
            case 5:
                mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); break;
            case 15:
                mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;
            //Early Exit
            case 6:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: ftol is too small (OPTI setting tolfun). No further improvement in the sum of squares is possible ***\n"); break;
            case 7:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: xtol is too small (OPTI hard-codes 1e-7). No further improvement in the approximate solution x is possible ***\n"); break;
            case 8:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: gtol is too small (OPTI hard-codes 1e-7). Fvec is orthogonal to the columns of the jacobian to machine precision ***\n"); break;
            case -1:
                mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
        }

        if(*exitflag==1)
            mexPrintf("\n Final SSE: %12.5g\n In %3.0f function evaluations\n",*fval,*iter);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Free Memory
    if(diag) mxFree(diag);
    if(fvec) mxFree(fvec);
    if(fjac) mxFree(fjac);
    if(ipvt) mxFree(ipvt);
    if(qtf)  mxFree(qtf);
    if(mptr) mxFree(mptr);
}

//LM_DIF Callback
static void func(int *m, int *n, double *x, double *fvec, int *iflag)
{
    bool havrnorm = false, stop = false;
    int i, stat;
    double *fval, rnorm = 0, evaltime;   
    
    fun.plhs[0] = NULL;
    memcpy(mxGetPr(fun.prhs[fun.xrhs]), x, *n * sizeof(double));
    
    stat = mexCallMATLAB(1, fun.plhs, fun.nrhs, fun.prhs, fun.f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");
    
    //Get Objective
    fval = mxGetPr(fun.plhs[0]);
    //Assign Objective + subtract ydata
    for(i=0;i<*m;i++)
        fvec[i] = fval[i]-ydata[i];
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Iteration Printing
    if(printLevel > 1) {
        if(*iflag==1) {//feval             
            //Calculate Residual Norm
            rnorm = 0; havrnorm = true;
            for(i=0;i<*m;i++)
                rnorm += (fval[i]-ydata[i])*(fval[i]-ydata[i]);

            if(citer == 1 || !(citer%10))
                mexPrintf(" feval       time            sse\n");

            mexPrintf("%5d       %5.2f     %12.5g\n",citer,evaltime,rnorm);
        }
        mexEvalString("drawnow;"); //flush draw buffer
    }
    
    //Iteration Callback
    if(iterF.enabled)
    {
        //Calculate sse if we don't have it
        if(!havrnorm)
            for(i=0;i<*m;i++)
                rnorm += (fval[i]-ydata[i])*(fval[i]-ydata[i]);

        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), &citer, sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), &rnorm, sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), x, *n * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
    }
    
    // Clean up Ptr
    mxDestroyArray(fun.plhs[0]);

    //Update feval counter
    citer++;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting LM_DIF...\n\n");
        *iflag = -1; //terminate
    }
    //Check for iterfun terminate
    if (stop) {
        mexPrintf("\nIterFun called Stop. Exiting LM_DER...\n\n");
        *iflag = -1; //terminate
    }
    //Check for maxtime expiry    
    if(evaltime > maxtime)
    {
        mexPrintf("\nMaximum Solver Time Exceeded. Exiting LM_DER...\n\n");
        *iflag = -1; //terminate
    }
}


//LM_DER Callback
static void jacfunc(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag)
{
    bool havrnorm = false, stop = false;
    int i, stat;
    double *fval, rnorm=0, evaltime; 
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    //Function Eval
    if(*iflag == 1) {   
        fun.plhs[0] = NULL;
        memcpy(mxGetPr(fun.prhs[fun.xrhs]), x, *n * sizeof(double));

        stat = mexCallMATLAB(1, fun.plhs, fun.nrhs, fun.prhs, fun.f);
        if(stat)
          mexErrMsgTxt("Error calling Objective Function!");
        //Get Objective
        fval = mxGetPr(fun.plhs[0]);
        //Assign Objective + subtract ydata
        for(i=0;i<*m;i++)
            fvec[i] = fval[i]-ydata[i];

        //Iteration Printing
        if(printLevel > 1) {      
            //Calculate Residual Norm
            rnorm = 0; havrnorm = true;
            for(i=0;i<*m;i++)
                rnorm += (fval[i]-ydata[i])*(fval[i]-ydata[i]);

            if(citer == 1 || !(citer%10))
                mexPrintf(" feval      time[s]          sse\n");

            mexPrintf("%5d       %5.2f     %12.5g\n",citer,evaltime,rnorm);
            mexEvalString("drawnow;"); //flush draw buffer
        }
        
        //Iteration Callback
        if(iterF.enabled)
        {
            //Calculate sse if we don't have it
            if(!havrnorm)
                for(i=0;i<*m;i++)
                    rnorm += (fval[i]-ydata[i])*(fval[i]-ydata[i]);
            
            iterF.plhs[0] = NULL;
            memcpy(mxGetData(iterF.prhs[1]), &citer, sizeof(int));
            memcpy(mxGetPr(iterF.prhs[2]), &rnorm, sizeof(double));
            memcpy(mxGetPr(iterF.prhs[3]), x, *n * sizeof(double));
            stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
            if(stat)
                mexErrMsgTxt("Error calling Callback Function!");
            
            //Collect return argument
            stop = *(bool*)mxGetData(iterF.plhs[0]);
            // Clean up Ptr
            mxDestroyArray(iterF.plhs[0]);
        }
        // Clean up Ptr
        mxDestroyArray(fun.plhs[0]);
        
        //Update feval counter
        citer++;
    }
    //Jacobian Eval
    else if(*iflag == 2) {
        grad.plhs[0] = NULL;
        memcpy(mxGetPr(grad.prhs[grad.xrhs]), x, *n * sizeof(double));

        stat = mexCallMATLAB(1, grad.plhs, grad.nrhs, grad.prhs, grad.f);
        if(stat)
          mexErrMsgTxt("Error calling Gradient Function!");

        //Get Objective
        fval = mxGetPr(grad.plhs[0]);
        //Assign Gradient
        memcpy(fjac,fval,*m**n*sizeof(double));

        // Clean up Ptr
        mxDestroyArray(grad.plhs[0]);
    }    
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting LM_DER...\n\n");
        *iflag = -1; //terminate
    }
    //Check for iterfun terminate
    if (stop) {
        mexPrintf("\nIterFun called Stop. Exiting LM_DER...\n\n");
        *iflag = -1; //terminate
    }
    //Check for maxtime expiry    
    if(evaltime > maxtime)
    {
        mexPrintf("\nMaximum Solver Time Exceeded. Exiting LM_DER...\n\n");
        *iflag = -1; //terminate
    }
}

void checkInputs(const mxArray *prhs[], int nrhs)
{    
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 4 arguments to lmder!\n\nlmder(fun,grad,x0,ydata)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(prhs[0]))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsEmpty(prhs[1]) && (!mxIsFunctionHandle(prhs[1]) && !mxIsChar(prhs[1])))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsEmpty(prhs[2]))
        mexErrMsgTxt("x0 must be a real double column vector!");
    if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsEmpty(prhs[3]))
        mexErrMsgTxt("ydata must be a real double column vector!");

    //Check Options
    if(nrhs > 4) {
        if(!mxIsStruct(prhs[4]))
            mexErrMsgTxt("The specified options must be a structure!");
    }

}

double getStatus(int stat)
{
    switch((int)stat)
    {        
        case 1:         //stopped by ftol
        case 2:         //stopped by xtol
        case 3:         //stopped by ftol and xtol
        case 4:         //stopped by gtol
            return 1;
            break;
        case 5:         //itmax
        case 15:        //max time
            return 0;
            break;
        case 6:         //ftol too small
        case 7:         //xtol too small
        case 8:         //gtol too small
            return -1;
        case 0:         //entry error
            return -2;
        case -1:
            return -5;  //user exit
        default:
            return -3;        
    }
}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" LM_DER + LM_DIF: MINPACK Levenberg-Marquardt Nonlinear Least Squares\n");           
    PRINT_BUILD_INFO;
    mexPrintf("  - Released as part of the MINPACK project: http://www.netlib.org/minpack/disclaimer\n");
    mexPrintf("  - Source available from: http://www.netlib.org/minpack/\n\n");

    mexPrintf(" MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}