/* PSWARMMEX - A MATLAB MEX Interface to PSWARM
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2013
 * www.i2c2.aut.ac.nz
 */

#include "mkl.h"
#include "mex.h"
#include "pswarm.h"
#include <time.h>
#include "opti_util.h"

#define PSWARM_VERSION "1.5"

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
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs, int *lincon);
double getStatus(int stat, int iter, int feval);
void func(int n, int m, double *x, double *lb, double *ub, double *fx);
double iterfcn(int n, int s, int iter, int gbest, struct swarm *pop);
extern int PSwarm(int, void (*)(),double *, double *, int, double *, double *, double **, double *, double *);

//Global Variables (no user arg passing to PSwarm functions...)
user_function_data fun;
iter_fun_data iterF;
int printLevel, noFeval;
double maxtime;
clock_t start, end;
bool ctrlCExit;
//Structures common to PSwarm and MEX
struct Options opt;
struct Stats stats;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args         
    double *x0 = NULL, *lb = NULL, *ub = NULL, *A = NULL, *b = NULL;  
    
    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval;
    
    //Internal Vars
    double *llb, *lub;
    size_t ndec;   
    int i, exit_code;    
    
    //PSwarm Vars
    double *sol = NULL;
    int lincon = 0;

    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(PSWARM_VERSION);
        return;
    }
    
    //Check user inputs
    checkInputs(prhs,nrhs,&lincon);
    
    //Set Defaults    
    printLevel = 0;
    maxtime = 1000;
    noFeval = 0;
    ctrlCExit = false;
    iterF.enabled = false;
    //Reset Stats
    stats.solveriters = 0;
    stats.objfunctions = 0;
    stats.pollsteps = 0;
    stats.sucpollsteps = 0;
    //Set PSwarm Defaults
    opt.s = 42;
    opt.mu = 0.5; opt.nu = 0.5;
    opt.maxvfactor = 0.5; opt.maxiter = 1500;
    opt.maxf = 10000;       
    opt.iweight = 0.9; opt.fweight = 0.4;
    opt.n2grd = 0.5;
    opt.blim = 10;
    opt.tol = 1e-5;
    opt.delta = Inf; opt.fdelta = 5.0; opt.idelta = 2.0; opt.ddelta = 0.5;
    opt.pollbasis = 0; opt.EpsilonActive = 0.1;
    opt.IPrint = 10;    
    opt.vectorized = 0;     //important, otherwise will pass whole swarm
    opt.outfcn = &iterfcn;  //iteration + ctrl c

    //Get Sizes
    ndec = mxGetNumberOfElements(prhs[1]);
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

    //Get x0
    x0 = mxGetPr(prhs[1]);   
    
    //Get Bounds
    //LB
    if(!mxIsEmpty(prhs[2])){
        llb = mxGetPr(prhs[2]);
        lb = mxCalloc(ndec,sizeof(double));
        memcpy(lb,llb,ndec*sizeof(double));
        for(i=0;i<ndec;i++) {
            if(mxIsInf(lb[i]))
                lb[i] = -1e19;
        }
    }
    else {
        lb = mxCalloc(ndec,sizeof(double));
        for(i=0;i<ndec;i++)
            lb[i] = -1e19;
    }
    //UB
    if(nrhs > 3 && !mxIsEmpty(prhs[3])){
        lub = mxGetPr(prhs[3]);
        ub = mxCalloc(ndec,sizeof(double));
        memcpy(ub,lub,ndec*sizeof(double));
        for(i=0;i<ndec;i++) {
            if(mxIsInf(ub[i]))
                ub[i] = 1e19;
        }
    }
    else {
        ub = mxCalloc(ndec,sizeof(double));
        for(i=0;i<ndec;i++)
            ub[i] = 1e19;
    }
    
    //Get Linear Inequality Constraints
    if(nrhs > 4) {
        if(!mxIsEmpty(prhs[4]) && !mxIsEmpty(prhs[5])) {
            A = mxGetPr(prhs[4]);
            b = mxGetPr(prhs[5]);
        }
    }
    
    //Get Options if specified
    if(nrhs > 6) {
        if(mxGetField(prhs[6],0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(prhs[6],0,"display"));
        if(mxGetField(prhs[6],0,"maxiter"))
            opt.maxiter = (int)*mxGetPr(mxGetField(prhs[6],0,"maxiter"));
        if(mxGetField(prhs[6],0,"maxtime"))
            maxtime = *mxGetPr(mxGetField(prhs[6],0,"maxtime"));
        if(mxGetField(prhs[6],0,"tolfun"))
            opt.tol = *mxGetPr(mxGetField(prhs[6],0,"tolfun"));
        if(mxGetField(prhs[6],0,"maxfeval"))
            opt.maxf = (int)*mxGetPr(mxGetField(prhs[6],0,"maxfeval"));
        if(mxGetField(prhs[6],0,"swarm_size"))
            opt.s = (int)*mxGetPr(mxGetField(prhs[6],0,"swarm_size"));
        if(mxGetField(prhs[6],0,"vectorized"))
            opt.vectorized = (int)*mxGetPr(mxGetField(prhs[6],0,"vectorized"));
        if(mxGetField(prhs[6],0,"mu"))
            opt.mu = *mxGetPr(mxGetField(prhs[6],0,"mu"));
        if(mxGetField(prhs[6],0,"nu"))
            opt.nu = *mxGetPr(mxGetField(prhs[6],0,"nu"));
        if(mxGetField(prhs[6],0,"iweight"))
            opt.iweight = *mxGetPr(mxGetField(prhs[6],0,"iweight"));
        if(mxGetField(prhs[6],0,"fweight"))
            opt.fweight = *mxGetPr(mxGetField(prhs[6],0,"fweight"));
        if(mxGetField(prhs[6],0,"delta"))
            opt.delta = *mxGetPr(mxGetField(prhs[6],0,"delta"));
        if(mxGetField(prhs[6],0,"idelta"))
            opt.idelta = *mxGetPr(mxGetField(prhs[6],0,"idelta"));
        if(mxGetField(prhs[6],0,"ddelta"))
            opt.ddelta = *mxGetPr(mxGetField(prhs[6],0,"ddelta"));
        if(mxGetField(prhs[6],0,"iterfun") && !mxIsEmpty(mxGetField(prhs[6],0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(prhs[6],0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
    }       
    
    //If not vectorized, we can create x now, otherwise must be done in callback
    if(!opt.vectorized)
        fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL);
    
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

    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is PSwarm v1.5\n");
            
        mexPrintf(" Authors: A.I.F Vaz and L.N. Vicente\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Linear Constraints:     %4d\n",lincon);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Run PSwarm
    start = clock();
    exit_code = PSwarm((int)ndec, &func, lb, ub, lincon, A, b, &sol, fval, x0);
    
    if(exit_code == 0) {
        //Copy Solution
        memcpy(x,sol,ndec*sizeof(double));
        free(sol);
        *iter = (double)stats.solveriters;
        *feval = (double)stats.objfunctions;
    }

    //Save Status & Iterations
    *exitflag = getStatus(exit_code,stats.solveriters,stats.objfunctions); 
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch((int)*exitflag)
        {
            //Success
            case 1:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Normal Exit ***\n"); break;
            //Error
            case -1:
                mexPrintf("\n *** ERROR: Abnormal Exit ***\n"); break;
            case 0:
                if(stats.solveriters >= (opt.maxiter-5))
                    mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); 
                else if(((double)(end-start))/CLOCKS_PER_SEC > maxtime)
                    mexPrintf("\n *** MAXIMUM TIME EXCEEDED ***\n");
                else
                    mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); 
                
                break;
            //Early Exit
            case -2:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Failed to allocate memory ***\n"); break;
            case -3:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Unable to initialize population - check constraints are feasible ***\n"); break;
            case -5:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: User Exit (Ctrl C) ***\n"); break;
        }

        if(*exitflag==1)
            mexPrintf("\n Final Objective: %12.5g\n In %3d iterations and\n   %4d function evaluations\n",*fval,stats.solveriters,stats.objfunctions);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Free Memory
    if(lb) mxFree(lb);
    if(ub) mxFree(ub);    
}
   
//PSwarm Callback
void func(int n, int m, double *x, double *lb, double *ub, double *fx)
{
    int stat;
    double *fval;   
    
    fun.plhs[0] = NULL;
    
    //If Vectorized create a new X, otherwise already done
    if(opt.vectorized)
        fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(n, m, mxREAL);
    //Copy x from PSwarm to MATLAB memory
    memcpy(mxGetPr(fun.prhs[fun.xrhs]), x, n * m * sizeof(double));

    stat = mexCallMATLAB(1, fun.plhs, fun.nrhs, fun.prhs, fun.f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");      
    
    if(opt.vectorized && mxGetNumberOfElements(fun.plhs[0]) != m)
        mexErrMsgTxt("Vectorization Error: The returned fval should be a column vector equal to the length of the swarm passed");
    
    //Get Objective and copy
    fval = mxGetPr(fun.plhs[0]);
    memcpy(fx,fval,m*sizeof(double));

    // Clean up Ptr
    mxDestroyArray(fun.plhs[0]);
    //Clean up x if vectorized
    if(opt.vectorized)
        mxDestroyArray(fun.prhs[fun.xrhs]);
    
    noFeval++;
}

double iterfcn(int n, int s, int iter, int gbest, struct swarm *pop)
{
    bool stop = false;
    int stat;
    double ret = 1.0;
    double evaltime;
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Iteration Printing
    if(printLevel > 1 && iter > 0) {      
        if(iter == 10 || !(iter%100))
            mexPrintf(" iter      feval     time[s]   leader     objective\n");
        
        mexPrintf("%5d     %5d     %5.2f      %5d  %12.5g\n",iter,noFeval,evaltime,gbest,pop->fy[gbest]);       
        mexEvalString("drawnow;"); //flush draw buffer        
    }
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting PSwarm...\n\n");
        ret = -1.0;
        ctrlCExit = true;
    }
    
    //Iteration Callback
    if(iterF.enabled && iter > 0)
    {
        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), &iter, sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), &pop->fy[gbest], sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), &pop->y[gbest*n], n * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        if(stop) 
        {   
            ctrlCExit = true;
            mexPrintf("\nIterFun Called Stop. Exiting PSwarm...\n\n");
            ret = -1;
        }
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
    }
    
    //Check for maxtime expiry    
    if(evaltime > maxtime)
    {
        mexPrintf("\nMaximum Solver Time Exceeded. Exiting PSwarm...\n\n");
        ret = -1.0;
    }   
    
    return ret;
}


void checkInputs(const mxArray *prhs[], int nrhs, int *lincon)
{    
    size_t Mx0;
    
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 4 arguments to pswarm!\n\npswarm(fun,x0,lb,ub)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(prhs[0]))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsEmpty(prhs[1]))
        mexErrMsgTxt("x0 must be a real double column vector!");

    //Get ndec
    Mx0 = mxGetNumberOfElements(prhs[2]);
    
    //Check Bounds
    if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgTxt("lb must be a real double column vector!");
    if(!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
        mexErrMsgTxt("ub must be a real double column vector!");
    //Check Empty
    if(mxIsEmpty(prhs[2]))
        mexErrMsgTxt("lb is empty! PSwarm requires finite bounds");
    if(mxIsEmpty(prhs[3]))
        mexErrMsgTxt("ub is empty! PSwarm requires finite bounds");
    //Check Sizes
    if(Mx0 != mxGetNumberOfElements(prhs[2]))
        mexErrMsgTxt("lb is not the same length as x0! Ensure they are both Column Vectors");
    if(Mx0 != mxGetNumberOfElements(prhs[3]))
        mexErrMsgTxt("ub is not the same length as x0! Ensure they are both Column Vectors");
    
    //Check Linear Inequality
    if(nrhs > 4) {
        if(nrhs == 5)
            mexErrMsgTxt("You must supply both A & b!");
        if(!mxIsEmpty(prhs[4]) && (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4])))
            mexErrMsgTxt("A must be a full, real double matrix!");
        if(!mxIsEmpty(prhs[5]) && (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])))
            mexErrMsgTxt("b must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(prhs[4])) {
            if(mxGetN(prhs[4]) != Mx0)
                mexErrMsgTxt("A is not the right size!");
            if(mxGetM(prhs[4]) != mxGetM(prhs[5]))
                mexErrMsgTxt("A & b sizes do not correspond");
            //Confirm no linear equations
            *lincon = (int)mxGetM(prhs[5]);
        }
    }
    
    //Check Options
    if(nrhs > 6) {
        if(!mxIsStruct(prhs[6]))
            mexErrMsgTxt("The specified options must be a structure!");
    }
}

double getStatus(int stat, int iter, int fevals)
{
    if(ctrlCExit)
        return -5;
    
    switch((int)stat)
    {        
        case 0:         //normal exit
            if(iter >= (opt.maxiter-5) || fevals >= (opt.maxf-5))
                return 0;
            else if(difftime(end,start) > maxtime)
                return 0;
            else
                return 1;
            break;
        case 1:         //abnormal exit
            return -1;
        case 2:         //memory error
            return -2;
        case 3:
            return -3;  //population error
        default:
            return -3;        
    }
}

//Print Solver Information
void printSolverInfo()
{    
    char vbuf[6]; getVSVer(vbuf);  
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" PSWARM: Particle Swarm and Pattern Based Optimization [v%s, Built %s, VS%s]\n",PSWARM_VERSION,__DATE__,vbuf);              
    mexPrintf("  - Released under the GNU Lesser General Public License: http://lpsolve.sourceforge.net/5.5/LGPL.htm\n");
    mexPrintf("  - Source available from: http://www.norg.uminho.pt/aivaz/pswarm/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}