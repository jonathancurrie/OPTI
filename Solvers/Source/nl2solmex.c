/* NL2SOLMEX - A MATLAB MEX Interface to NL2SOL
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2011-2013
 * www.i2c2.aut.ac.nz
 */

#include "mex.h"
#include "mkl.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

#define NL2SOL_VERSION "2.3"

//Argument Enumeration (in expected order of arguments)
enum {eFUN, eGRAD, eX0, eYDATA, eLB, eUB, eOPTS};
//PRHS Defines    
#define pFUN    prhs[eFUN]
#define pGRAD   prhs[eGRAD]
#define pX0     prhs[eX0]
#define pYDATA  prhs[eYDATA]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pOPTS   prhs[eOPTS]

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN], g[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     mxArray *prhs_g[MAXRHS];
     int xrhs, nrhs, xrhs_g, nrhs_g;
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

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
double getStatus(int stat);
static void CALCR(int *n, int *p, double *x, int *nf, double *r, int *uiparm, double *ydata, void *ufparm);
static void CALCJ(int *n, int *p, double *x, int *nf, double *j, int *uiparm, double *ydata, void *ufparm);

//NEW ROUTINES (v2.3)
// DN2G routine (implements NL2SOL v2.3)
extern void DN2G(int *n, int *p, double *x, 
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 void(*jac)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);
// DN2GB routine 
extern void DN2GB(int *n, int *p, double *x, double *b,
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 void(*jac)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

// DN2F routine (implements NL2SNO v2.3)
extern void DN2F(int *n, int *p, double *x, 
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);
// DN2FB routine 
extern void DN2FB(int *n, int *p, double *x, double *b,
                 void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                 int *iv, int *liv, int *lv, double *v, int *uiparm, double *urparm, void *ufparm);

// DIVSET routine to initialize settings
extern void DIVSET (int *alg, int *iv, int *liv, int *lv, double *v);
extern double D1MACH(int *val);

//OLD ROUTINES (v2.2)
/*// NL2SOL routine implemented in Fortran 77.
extern void NL2SOL (int *n, int *p, double *x, 
                    void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*), 
                    void(*jac)(int*,int*,double*,int*,double*,int*,double*,void*), 
                    int *iv, double *v, int *uiparm, double *urparm, void *ufparm);

// NL2SNO routine implemented in Fortran 77. 
extern void NL2SNO (int *n, int *p, double *x, 
                    void(*fun)(int*,int*,double*,int*,double*,int*,double*,void*),
                    int *iv, double *v, int *uiparm, double *urparm, void* ufparm);

// DFAULT routine to initialise settings
extern void DFAULT (int *iv, double *v);*/

//Iteration Callback Structure
iter_fun_data iterF;
//Timing Information
clock_t start;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    user_function_data fun;
    double *x0;          
    double *ydata;  
    double *lb, *ub;

    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval;
    
    //Internal Vars
    size_t ndec, ndat, i, j;   
    int havJac = 0;
    int printLevel = 0;
    double *bounds = NULL;
    
    //NL2SOL Vars
    int n;                  //len data
    int p;                  //len x    
    int *iv;                //intermediate work array
    double *v;              //intermediate work array
    int liv = 0, lv = 0;    //lengths of intermediate arrays
    int uiparm[2] = {1,0};  //user integer array [citer, printLevel]
    int one = 1, two = 2;
    //Defaults
    int maxfev = 1500;
    int maxiter = 1000;
    double frtol = 1e-7;
    double fatol = 1e-5;
    iterF.enabled = false;

    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(NL2SOL_VERSION);
            
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs);

    //Get Sizes
    ndec = mxGetNumberOfElements(pX0);
    ndat = mxGetNumberOfElements(pYDATA);
    //Get Objective Function Handle
    if (mxIsChar(pFUN)) {
        CHECK(mxGetString(pFUN, fun.f, FLEN) == 0,"error reading objective name string");
        fun.nrhs = 1;
        fun.xrhs = 0;
    } else {
        fun.prhs[0] = (mxArray*)pFUN;
        strcpy(fun.f, "feval");
        fun.nrhs = 2;
        fun.xrhs = 1;
    }
    fun.prhs[fun.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    //Check and Get Gradient Function Handle 
    if(!mxIsEmpty(pGRAD)) {  
        havJac = 1;
        if (mxIsChar(pGRAD)) {
            CHECK(mxGetString(pGRAD, fun.g, FLEN) == 0,"error reading gradient name string");
            fun.nrhs_g = 1;
            fun.xrhs_g = 0;
        } else {
            fun.prhs_g[0] = (mxArray*)pGRAD;
            strcpy(fun.g, "feval");
            fun.nrhs_g = 2;
            fun.xrhs_g = 1;
        }   
        fun.prhs_g[fun.xrhs_g] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    }    

    //Get x0 + data
    x0 = mxGetPr(pX0);
    ydata = mxGetPr(pYDATA);   
    
    //Get bounds if specified
    if(nrhs > eUB && !mxIsEmpty(pUB))
    {
        lb = mxGetPr(pLB);
        ub = mxGetPr(pUB);
        bounds = (double*)mxCalloc(2*ndec,sizeof(double));
        //Copy bounds into nlsol b(2,p) array
        for(i=0,j=0;i<2*ndec;i+=2,j++)
        {
            if(mxIsInf(lb[j]))
                bounds[i] = -D1MACH(&two);
            else
                bounds[i] = lb[j];
            if(mxIsInf(ub[j]))
                bounds[i+1] = D1MACH(&two);
            else
                bounds[i+1] = ub[j];
        }
    }
    
    //Get Options if specified
    if(nrhs > eOPTS) {
        if(mxGetField(pOPTS,0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display"));
        if(mxGetField(pOPTS,0,"maxfeval"))
            maxfev = (int)*mxGetPr(mxGetField(pOPTS,0,"maxfeval"));
        if(mxGetField(pOPTS,0,"maxiter"))
            maxiter = (int)*mxGetPr(mxGetField(pOPTS,0,"maxiter"));
        if(mxGetField(pOPTS,0,"tolrfun"))
            frtol = *mxGetPr(mxGetField(pOPTS,0,"tolrfun"));
        if(mxGetField(pOPTS,0,"tolafun"))
            fatol = *mxGetPr(mxGetField(pOPTS,0,"tolafun"));
        if(mxGetField(pOPTS,0,"iterfun") && !mxIsEmpty(mxGetField(pOPTS,0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(pOPTS,0,"iterfun");
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
    
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        if(havJac) {
            if(bounds)
                mexPrintf(" This is DN2GB (NL2SOL v%s)\n",NL2SOL_VERSION);
            else
                mexPrintf(" This is DN2G (NL2SOL v%s)\n",NL2SOL_VERSION);
        }
        else {
            if(bounds)
                mexPrintf(" This is DN2FB (NL2SNO v%s)\n",NL2SOL_VERSION);
            else
                mexPrintf(" This is DN2F (NL2SNO v%s)\n",NL2SOL_VERSION);
        }
            
        mexPrintf(" Authors: John Dennis, David Gay, Roy Welsch\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Data Points:            %4d\n",ndat);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Assign Arguments
    n = (int)ndat;
    p = (int)ndec;
    if(bounds)
        liv = 82 + 4*p;
    else
        liv = 82+p;
    if(bounds)
        lv = 105 + p*(n + 2*p + 21) + 2*n;
    else
        lv = 105 + p*(n + 2*p + 17) + 2*n;
    iv = (int*)mxCalloc(liv + 10,sizeof(int)); 
    v  = (double*)mxCalloc(lv + 10,sizeof(double));
    //Set Options
    //DFAULT(iv,v); //set defaults (OLD v2.2)
    DIVSET(&one,iv,&liv,&lv,v);
    iv[13] = iv[14] = 0;        //no covar
    iv[16] = maxfev;    //limit on fevals + gevals
    iv[17] = maxiter;   //max iter
    iv[18] = 0;     //no iteration printing
    iv[19] = 0;     //no default printing
    iv[20] = 0;     //no output unit printing
    iv[21] = 0;     //no x printing
    iv[22] = 0;     //no summary printing
    iv[23] = 0;     //no initial printing
    v[30] = fatol;   
    v[31] = frtol;
    //MEX Options
    uiparm[1] = printLevel;
    
    //Start timer
    start = clock();
    
    //Call Algorithm based on Derivatives    
    if(havJac)
    {
        //NL2SOL(&n,&p,x,CALCR,CALCJ,iv,v,uiparm,ydata,&fun); 
        if(bounds)
            DN2GB(&n,&p,x,bounds,CALCR,CALCJ,iv,&liv,&lv,v,uiparm,ydata,&fun);
        else
            DN2G(&n,&p,x,CALCR,CALCJ,iv,&liv,&lv,v,uiparm,ydata,&fun);
    }
    else 
    {
        //#ifdef WIN64
        //    mexErrMsgTxt("Currently NL2SNO crashes when compiled on 64bit systems, try NL2SOL instead (supply a gradient)");    
        //#else
            //NL2SNO(&n,&p,x,CALCR,iv,v,uiparm,ydata,&fun);
        if(bounds)
            DN2FB(&n,&p,x,bounds,CALCR,iv,&liv,&lv,v,uiparm,ydata,&fun);
        else
            DN2F(&n,&p,x,CALCR,iv,&liv,&lv,v,uiparm,ydata,&fun);
        //#endif
    }
    
    //Final Rnorm
    *fval = 2*v[9]; //nl2sol uses 1/2 sum(resid^2)
    //Save Status & Iterations
    *exitflag = getStatus(iv[0]);
    *iter = (double)iv[30];
    *feval = (double)uiparm[0];
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch((int)iv[0])
        {
            //Success
            case 3:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** x-convergence ***\n"); break;
            case 4:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** relative function convergence ***\n"); break;
            case 5:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** both x and relative function convergence ***\n"); break;
            case 6:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** absolute function convergence ***\n"); break;
            //Error
            case 9:
                mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); break;
            case 10:
                mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;                       
            case 13:
                mexPrintf("\n *** ERROR: f(x) cannot be computed at the initial x ***\n"); break;
            case 14:
                mexPrintf("\n *** ERROR: bad parameters passed to asses ***\n"); break;    
            case 15:
                mexPrintf("\n *** ERROR: the jacobian cannot be computed at x ***\n"); break;       
            case 16:
                mexPrintf("\n *** ERROR: internal parameter n or p out of range ***\n"); break;
            case 17:
                mexPrintf("\n *** ERROR: restart attempted with n or p changed ***\n"); break;    
            case 50:
                mexPrintf("\n *** ERROR: internal parameter iv(1) is out of range ***\n"); break;
            //Early Exit
            case 7:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: singular convergence. the hessian near the current iterate appears to be singular or nearly so. ***\n"); break;
            case 8:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: false convergence. the iterates appear to be converging to a noncritical point. this may mean that the convergence tolerances are too small ***\n"); break;
            case 11:
                mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
            //Other Error
            default:
                mexPrintf("\n *** ERROR: internal error code %d ***\n",(int)iv[0]); break;
        }
        
        if(*exitflag==1)
            mexPrintf("\n Final SSE: %12.5g\n In %3.0f iterations\n",*fval,*iter);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Free Memory
    if(iv) mxFree(iv); iv = NULL;
    if(v) mxFree(v); v = NULL;
    if(bounds) mxFree(bounds); bounds = NULL;
}

//CALCR
static void CALCR(int *n, int *p, double *x, int *nf, double *r, int *uiparm, double *ydata, void *ufparm)
{
    bool havrnorm = false, stop = false;
    int i, stat;
    double *fval, rnorm = 0;  
    clock_t end;
    double evaltime;
    
    //Get User Data
    user_function_data *fun = (user_function_data*)ufparm;
    
    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, *p * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");
    
    //Get Objective
    fval = mxGetPr(fun->plhs[0]);

    //Assign Objective + subtract ydata
    for(i=0;i<*n;i++) {
        r[i] = fval[i]-ydata[i];

        //Check for out of bounds, will try a smaller step size
        if(mxIsInf(r[i]) || mxIsNaN(r[i]))
            *nf = 0;
    }
    
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
    
    //Iteration Printing
    if(uiparm[1] > 1) {   
        //Get Execution Time
        end = clock();
        evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
        //Calculate Residual Norm
        rnorm = 0; havrnorm = true;
        for(i=0;i<*n;i++)
            rnorm += r[i]*r[i];
        
        if(uiparm[0] == 1 || !(uiparm[0]%10))
            mexPrintf(" feval      time[s]          sse\n");

        mexPrintf("%5d       %5.2f     %12.5g\n",uiparm[0],evaltime,rnorm);
        mexEvalString("drawnow;"); //flush draw buffer
    }
    
    //Iteration Callback
    if(iterF.enabled)
    {
        //Calculate sse if we don't have it
        if(!havrnorm)
            for(i=0;i<*n;i++)
                rnorm += r[i]*r[i];;

        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), uiparm, sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), &rnorm, sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), x, *p * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
        
        //Warn user stop not implemented
        if(stop)
            mexWarnMsgTxt("NL2SOL does not implement the stop feature of iterfun");
    }
    
    uiparm[0]++;
}

//CALCJ
static void CALCJ(int *n, int *p, double *x, int *nf, double *j, int *uiparm, double *ydata, void *ufparm)
{
    int stat;
    double *grad;   
    
    user_function_data *fun = (user_function_data *)ufparm;

    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs_g[fun->xrhs_g]), x, *p * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs_g, fun->prhs_g, fun->g);
    if(stat)
      mexErrMsgTxt("Error calling Gradient Function!");
    
    //Get Gradient
    grad = mxGetPr(fun->plhs[0]);
    //Assign Gradient
    memcpy(j,grad,*n**p*sizeof(double));
    
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
}

void checkInputs(const mxArray *prhs[], int nrhs)
{    
    size_t ndec;
    
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 4 arguments to nl2sol!\n\nnl2sol(fun,grad,x0,ydata)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(pFUN))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsEmpty(pGRAD) && (!mxIsFunctionHandle(pGRAD) && !mxIsChar(pGRAD)))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");
    if(!mxIsDouble(pYDATA) || mxIsComplex(pYDATA) || mxIsEmpty(pYDATA))
        mexErrMsgTxt("ydata must be a real double column vector!");

    ndec = mxGetNumberOfElements(pX0);
    
    //Check for lb + ub
    if(nrhs == eUB)
    	mexErrMsgTxt("You must supply both lb and ub for bounded problems!");
    
    //Check Bounds
    if(nrhs > eUB)
    {
        if(mxIsEmpty(pLB) ^ mxIsEmpty(pUB))
            mexErrMsgTxt("You must supply both lb and ub for bounded problems");        
        if(!mxIsEmpty(pLB)) {
            if(!mxIsDouble(pLB) || mxIsComplex(pLB))
                mexErrMsgTxt("lb must be a real double column vector!");
            if(mxGetNumberOfElements(pLB) != ndec)
                mexErrMsgTxt("The number of elements in lb does not match x0");
        }
        if(!mxIsEmpty(pUB)) {
            if(!mxIsDouble(pUB) || mxIsComplex(pUB))
                mexErrMsgTxt("ub must be a real double column vector!");
            if(mxGetNumberOfElements(pUB) != ndec)
                mexErrMsgTxt("The number of elements in ub does not match x0");        
        }
    }
    
    //Check Options
    if(nrhs > eOPTS) {
        if(!mxIsStruct(pOPTS))
            mexErrMsgTxt("The specified options must be a structure!");
    }

}

double getStatus(int stat)
{
    switch((int)stat)
    {     
        case 3:         //stopped by xtol
        case 4:         //stopped by ftol
        case 5:         //stopped by ftol and xtol
        case 6:         //stopped by abs ftol
            return 1;
            break;
        case 9:         //feval max
        case 10:        //itmax
            return 0;
            break;
        case 8:         //tol too small
            return -1;
        case 13:        //bad initial f(x)
        case 14:        //bad parameters
        case 15:        //bad g(x)
        case 16:        //n or p out of range
        case 17:        //restart attempted (?)
        case 18:        //iv out of range
        case 19:        //v out of range
        case 50:        //iv[0] out of range
        case 87:        //v problem
            return -2;
        case 11:
            return -5;  //user exit
        default:
            return -3;        
    }
}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" NL2SOL: Adaptive Nonlinear Least Squares [v%s, Built %s]\n",NL2SOL_VERSION,__DATE__);              
    mexPrintf("  - Source available from: http://people.sc.fsu.edu/~jburkardt/f_src/nl2sol/nl2sol.html\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}