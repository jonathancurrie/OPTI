/* M1QN3MEX - A MATLAB MEX Interface to M1QN3
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.i2c2.aut.ac.nz
 */

#include "mex.h"
#include "mkl.h"
#include <stdio.h>
#include <string.h>
#include <time.h>

#define M1QN3_VERSION "3.3"

//Argument Enumeration (in expected order of arguments)
enum {eFUN, eGRAD, eX0, eOPTS};
//PRHS Defines    
#define pFUN    prhs[eFUN]
#define pGRAD   prhs[eGRAD]
#define pX0     prhs[eX0]
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

//Ctrl-C Detection
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);

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
static void SIMUL(int *indic, int *n, double *x, double *f, double *g, int *izs, float *rzs, double *dzs);

//M1QN3 Routine
extern void M1QN3(void(*fun)(int*,int*,double*,double*,double*,int*,float*,double*),
                  void(*prosca)(int*,double*,double*,double*,int*,float*,double*),
                  void(*ctonb)(int*,double*,double*,int*,float*,double*),
                  void(*ctcab)(int*,double*,double*,int*,float*,double*),  
                  int *n, double *x, double *f, double *g, double *dxmin, double *df1, double *epsg, 
                  char *normtype, int *impres, int *io, int *imode, int *omode, int *niter, int *nsim,
                  int *iz, double *dz, int *ndz, int *reverse, int *indec, int *izs, float *rzs, double *dzs);

//DEFAULT Routines
extern void EUCLID(int*,double*,double*,double*,int*,float*,double*);
extern void CTONBE(int*,double*,double*,int*,float*,double*);
extern void CTCABE(int*,double*,double*,int*,float*,double*);

//User Function Structure
user_function_data fun;
//Iteration Callback Structure
iter_fun_data iterF;
//Max Time data
double maxtime;
clock_t start, end;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    double *x0;          

    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval;
    
    //Internal Vars
    size_t ndec;  
    int printLevel = 0;
    
    //M1QN3 Vars
    int m = 5, n, indic, ndz;
    int uiparm[2] = {1,0};  //user integer array [citer, printLevel]
    double *g, dxmin = 1e-8, df1, epsg = 1e-6, *dz;
    char *normtype = "dfn";
    int impres = 0, io = 0, omode = 0, reverse = 0;
    int imode[3] = {0,0,0}; //DIS, cold start, no SIMUL with indic = 1
    int iz[5];
    float *rzs = NULL; double *dzs = NULL; //dummy args
    //Defaults
    int maxfev = 1500;
    int maxiter = 1000;
    maxtime = 1000;
    iterF.enabled = false;

    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(M1QN3_VERSION);
            
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs);

    //Get Sizes
    ndec = mxGetNumberOfElements(pX0);
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
    //Get Gradient Function Handle 
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

    //Get x0
    x0 = mxGetPr(pX0);
    
    //Get Options if specified
    if(nrhs > eOPTS) {
        if(mxGetField(pOPTS,0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display"));
        if(mxGetField(pOPTS,0,"maxfeval"))
            maxfev = (int)*mxGetPr(mxGetField(pOPTS,0,"maxfeval"));
        if(mxGetField(pOPTS,0,"maxiter"))
            maxiter = (int)*mxGetPr(mxGetField(pOPTS,0,"maxiter"));
        if(mxGetField(pOPTS,0,"maxtime"))
            maxtime = *mxGetPr(mxGetField(pOPTS,0,"maxtime"));
        if(mxGetField(pOPTS,0,"tolafun"))
            epsg = *mxGetPr(mxGetField(pOPTS,0,"tolafun")); //not function tolerance (gradient)
        if(mxGetField(pOPTS,0,"nupdates"))
            m = (int)*mxGetPr(mxGetField(pOPTS,0,"nupdates")); //number of l-bfgs updates
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
        mexPrintf(" This is M1QN3 v%s\n",M1QN3_VERSION);  
        mexPrintf(" Authors: Jean Charles Gilbert, Claude Lemarechal, INRIA\n MEX Interface J. Currie 2012\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Assign Arguments
    n = (int)ndec;
    indic = 4;
    g = (double*)mxCalloc(n,sizeof(double)); //allocate memory for gradient
    ndz = 4*n + m*(2*n + 1);
    dz = (double*)mxCalloc(ndz,sizeof(double));
    
    //Start timer
    start = clock();
    
    //Initialization Call
    SIMUL(&indic, &n, x, fval, g, uiparm, NULL, NULL);    
    //Set df1 (initial estiamte of f reduction)
    df1 = *fval;
    
    //MEX Options
    uiparm[0] = 1;
    uiparm[1] = printLevel;
    
    //Call Algorithm
    M1QN3(SIMUL,EUCLID,CTONBE,CTCABE,&n,x,fval,g,&dxmin,&df1,&epsg,normtype,
          &impres,&io,imode,&omode,&maxiter,&maxfev,iz,dz,&ndz,&reverse,&indic,
          uiparm,rzs,dzs);

    //Save Status & Iterations
    *exitflag = (double)omode;
    *iter = maxiter;
    *feval = maxfev;
    
    //Check if maxtime exceeded
    if(((double)(end-start))/CLOCKS_PER_SEC > maxtime)
        *exitflag = 8;
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch((int)*exitflag)
        {
            //Success
            case 1:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** gradient convergence |gk|/|g1| < epsg ***\n"); break;
            //Error
            case 5:
                mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); break;
            case 4:
                mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;                       
            case 8:
                mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;  
            case 2:
                mexPrintf("\n *** ERROR: one of the input arguments is not well initialized ***\n"); break;
            case 3:
                mexPrintf("\n *** ERROR: the line-search is blocked on tmax = 10^20 ***\n"); break;    
            case 6:
                mexPrintf("\n *** ERROR: stop dxmin during the line-search ***\n"); break;       
            case 7:
                mexPrintf("\n *** ERROR: either (g,d) is nonnegative or (y,s) is nonpositive ***\n"); break;
            //Early Exit
            case 0:
                mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
            //Other Error
            default:
                mexPrintf("\n *** ERROR: internal error code %d ***\n",omode); break;
        }
        
        if(*exitflag==1)
            mexPrintf("\n Final fval: %12.5g\n In %3.0f iterations\n",*fval,*iter);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Free Memory
    mxFree(g);
    mxFree(dz);
}

//MATLAB Callback
static void SIMUL(int *indic, int *n, double *x, double *f, double *g, int *izs, float *rzs, double *dzs)
{
    int stat;
    double *grad;   
    bool stop = false;
    double evaltime;

    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting M1QN3...\n\n");
        *indic = 0; //terminate
        return;
    }
    
    //Check for maxtime expiry    
    if(evaltime > maxtime)
    {
        mexPrintf("\nMaximum Solver Time Exceeded. Exiting M1QN3...\n\n");
        *indic = 0; //terminate
        return;
    }
    
    //Only compute f and g if requested
    if(*indic == 4)
    {    
        fun.plhs[0] = NULL;
        memcpy(mxGetPr(fun.prhs[fun.xrhs]), x, *n * sizeof(double));

        stat = mexCallMATLAB(1, fun.plhs, fun.nrhs, fun.prhs, fun.f);
        if(stat)
          mexErrMsgTxt("Error calling Objective Function!");

        //Get Objective
        *f = *mxGetPr(fun.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(fun.plhs[0]);
        
        //Check for inf, nan
        if(mxIsInf(*f) || mxIsNaN(*f))
            *indic = -1; //indicate smaller step size

        //Get Gradient
        fun.plhs[0] = NULL;
        memcpy(mxGetPr(fun.prhs_g[fun.xrhs_g]), x, *n * sizeof(double));

        stat = mexCallMATLAB(1, fun.plhs, fun.nrhs_g, fun.prhs_g, fun.g);
        if(stat)
          mexErrMsgTxt("Error calling Gradient Function!");

        //Get Gradient
        grad = mxGetPr(fun.plhs[0]);
        //Assign Gradient
        memcpy(g,grad,*n*sizeof(double));

        // Clean up Ptr
        mxDestroyArray(fun.plhs[0]);

        //Iteration Printing
        if(izs[1] > 1) {               
            if(izs[0] == 1 || !(izs[0]%10))
                mexPrintf(" feval       time           fval\n");

            mexPrintf("%5d       %5.2f    %12.5g\n",izs[0],evaltime,*f);
            mexEvalString("drawnow;"); //flush draw buffer
        }

        //Iteration Callback
        if(iterF.enabled)
        {
            iterF.plhs[0] = NULL;
            memcpy(mxGetData(iterF.prhs[1]), izs, sizeof(int));
            memcpy(mxGetPr(iterF.prhs[2]), f, sizeof(double));
            memcpy(mxGetPr(iterF.prhs[3]), x, *n * sizeof(double));
            stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
            if(stat)
                mexErrMsgTxt("Error calling Callback Function!");

            //Collect return argument
            stop = *(bool*)mxGetData(iterF.plhs[0]);
            // Clean up Ptr
            mxDestroyArray(iterF.plhs[0]);

            if(stop)
                *indic = 0; //force exit
        }
    
        //Increment feval counter
        izs[0]++;
    }
}

void checkInputs(const mxArray *prhs[], int nrhs)
{    
    if(nrhs < 3)
        mexErrMsgTxt("You must supply at least 3 arguments to m1qn3!\n\nm1qn3(fun,grad,x0)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(prhs[0]) && !mxIsChar(pFUN))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsFunctionHandle(pGRAD) && !mxIsChar(pGRAD))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");

    //Check Options
    if(nrhs > eOPTS) {
        if(!mxIsStruct(pOPTS))
            mexErrMsgTxt("The specified options must be a structure!");
    }

}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" M1QN3: Large-Scale Unconstrained L-BFGS Minimization [Built %s]\n",__DATE__);              
    mexPrintf("  - Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html\n");
    mexPrintf("  - Source available from: https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}
