/* MKLTRNLSMEX - A MATLAB MEX Interface to MKL Trust Region NLS
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.inverseproblem.co.nz
 */

#include "mex.h"
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "opti_build_utils.h"

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
void checkInputs(const mxArray *prhs[], int nrhs, int *con);
double getStatus(int stat);
int calcF(MKL_INT *ndec, MKL_INT *ndat, double *x, double *f, double *ydata, void *userdata, int *uiparm);
void calcG(MKL_INT *ndec, MKL_INT *ndat, double *x, double *g, void *userdata);
void calcNDG(MKL_INT *m, MKL_INT *n, double *x, double *f, void *data);

//Iteration callback
iter_fun_data iterF;
//Max Time data
double maxtime;
clock_t start, end;

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Input Args
    user_function_data fun;
    user_function_data grad;
    double *x0;          
    double *ydata = NULL, *lb = NULL, *ub = NULL;    

    //Outputs Args
    double *x, *fval, *exitflag, *iter, *feval, *geval;
    
    //Internal Vars  
    int havJac = 0, jacEval = 0;
    int uiparm[2] = {1,0};  //user integer array [citer, printLevel]
    int con = 0;
    int exit = 0;
    double *llb, *lub;
    
    //MKL Vars
    MKL_INT ndec, ndat;
    MKL_INT RCI_Request;
    MKL_INT successful, res;
    MKL_INT i;
    MKL_INT miter, st_cr;
    double *fvec = NULL, *fjac = NULL;
    double r1, r2;
    _TRNSP_HANDLE_t handle;
    
    //Defaults
    double eps[6] = {1e-6,1e-6,1e-6,1e-6,1e-6,1e-10};    //precisions for stop criteria
    double rs = 0.1;                                    //initial step bound
    MKL_INT iter1 = 1000, iter2 = 200;
    double tol = 1e-6;                                  //djacobi tol
    maxtime = 1000;
    iterF.enabled = false;

    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            char verstr[128];
            sprintf(verstr,"%d.%d R%d",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
            plhs[0] = mxCreateString(verstr);
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs,&con);

    //Get Sizes
    ndec = (MKL_INT)mxGetNumberOfElements(prhs[2]);
    ndat = (MKL_INT)mxGetNumberOfElements(prhs[3]);
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
    
    //Get Bounds
    if(con) {
        //LB
        if(!mxIsEmpty(prhs[4])){
            llb = mxGetPr(prhs[4]);
            lb = mxCalloc(ndec,sizeof(double));
            memcpy(lb,llb,ndec*sizeof(double));
            for(i=0;i<ndec;i++) {
                if(mxIsInf(lb[i]))
                    lb[i] = -1e30;
            }
        }
        else {
            lb = mxCalloc(ndec,sizeof(double));
            for(i=0;i<ndec;i++)
                lb[i] = -1e30;
        }
        //UB
        if(nrhs > 5 && !mxIsEmpty(prhs[5])){
            lub = mxGetPr(prhs[5]);
            ub = mxCalloc(ndec,sizeof(double));
            memcpy(ub,lub,ndec*sizeof(double));
            for(i=0;i<ndec;i++) {
                if(mxIsInf(ub[i]))
                    ub[i] = 1e30;
            }
        }
        else {
            ub = mxCalloc(ndec,sizeof(double));
            for(i=0;i<ndec;i++)
                ub[i] = 1e30;
        }
    }         
    
    //Get Options if specified
    if(nrhs > 6) {
        if(mxGetField(prhs[6],0,"display") && !mxIsEmpty(mxGetField(prhs[6],0,"display")))
            uiparm[1] = (int)*mxGetPr(mxGetField(prhs[6],0,"display"));
        if(mxGetField(prhs[6],0,"maxiter") && !mxIsEmpty(mxGetField(prhs[6],0,"maxiter")))
            iter1 = (int)*mxGetPr(mxGetField(prhs[6],0,"maxiter"));
        if(mxGetField(prhs[6],0,"maxtime") && !mxIsEmpty(mxGetField(prhs[6],0,"maxtime")))
            maxtime = *mxGetPr(mxGetField(prhs[6],0,"maxtime"));
        if(mxGetField(prhs[6],0,"iterfun") && !mxIsEmpty(mxGetField(prhs[6],0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(prhs[6],0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
        if(mxGetField(prhs[6],0,"tolTR") && !mxIsEmpty(mxGetField(prhs[6],0,"tolTR")))
            eps[0] = *mxGetPr(mxGetField(prhs[6],0,"tolTR"));
        if(mxGetField(prhs[6],0,"tolafun") && !mxIsEmpty(mxGetField(prhs[6],0,"tolafun")))
            eps[1] = *mxGetPr(mxGetField(prhs[6],0,"tolafun"));
        if(mxGetField(prhs[6],0,"tolSNGLR") && !mxIsEmpty(mxGetField(prhs[6],0,"tolSNGLR")))
            eps[2] = *mxGetPr(mxGetField(prhs[6],0,"tolSNGLR"));
        if(mxGetField(prhs[6],0,"tolSTEP") && !mxIsEmpty(mxGetField(prhs[6],0,"tolSTEP")))
            eps[3] = *mxGetPr(mxGetField(prhs[6],0,"tolSTEP"));
        if(mxGetField(prhs[6],0,"tolrfun") && !mxIsEmpty(mxGetField(prhs[6],0,"tolrfun")))
            eps[4] = *mxGetPr(mxGetField(prhs[6],0,"tolrfun"));
        if(mxGetField(prhs[6],0,"tolSTEPprec") && !mxIsEmpty(mxGetField(prhs[6],0,"tolSTEPprec")))
            eps[5] = *mxGetPr(mxGetField(prhs[6],0,"tolSTEPprec"));
        if(mxGetField(prhs[6],0,"trialStepSize") && !mxIsEmpty(mxGetField(prhs[6],0,"trialStepSize")))
            rs = *mxGetPr(mxGetField(prhs[6],0,"trialStepSize"));
        
        CheckOptiVersion(prhs[6]);
    }       
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,1, mxREAL);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    iter = mxGetPr(plhs[3]);
    feval = mxGetPr(plhs[4]);
    geval = mxGetPr(plhs[5]);
    
    //Copy initial guess to x
    memcpy(x,x0,ndec*sizeof(double));
    
    //Print Header
    if(uiparm[1]) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is the Intel MKL Trust Region NLS Solver\n");
        mexPrintf(" Author: Intel Math Kernel Library\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Data Points:            %4d\n",ndat);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Allocate memory and set initial values
    fvec = (double*) malloc (sizeof (double)*ndat);
	fjac = (double*) malloc (sizeof (double)*ndat*ndec);
	for (i = 0; i < ndat; i++)
		fvec [i] = 0.0;
	for (i = 0; i < ndat*ndec; i++)
		fjac [i] = 0.0;
    
    //Start timer
    start = clock();
    
    //Initialize Solver
    if(con)
        res = dtrnlspbc_init(&handle,&ndec,&ndat,x,lb,ub,eps,&iter1,&iter2,&rs);
    else
        res = dtrnlsp_init(&handle,&ndec,&ndat,x,eps,&iter1,&iter2,&rs);
        
    if(res != TR_SUCCESS)
    {
        *exitflag = -4;
        MKL_Free_Buffers();
        return;
    }
    //Set Initial RCI
    RCI_Request = 0;
    successful = 0;
    
    //Main RCI Cycle
    while(successful == 0)
    {
        //Call the solver
        if(con)
            res = dtrnlspbc_solve(&handle,fvec,fjac,&RCI_Request);
        else
            res = dtrnlsp_solve(&handle,fvec,fjac,&RCI_Request);
            
        //mexPrintf("solve res: %d, RCI: %d\n",res,RCI_Request);
        
        if(res != TR_SUCCESS)
        {
            *exitflag = -5;
            MKL_Free_Buffers();
            return;
        }
        //Determine we do next step
        if (RCI_Request == -1 || RCI_Request == -2 || RCI_Request == -3 || RCI_Request == -4 || RCI_Request == -5 || RCI_Request == -6)
			successful = 1; //exit

        if(RCI_Request == 1) { //recalculate function value
            exit = calcF(&ndec,&ndat,x,fvec,ydata,&fun,uiparm);
            if(exit < 1) { //user requested exit
               //Calculate final residual
                *fval = 0;
                for(i = 0; i < ndat; i++)
                    *fval += fvec[i]*fvec[i]; 
                *exitflag = -5;
                MKL_Free_Buffers();
                return;
            }
        }
        
        if(RCI_Request == 2) //compute gradient
        {
            jacEval++;
            if(havJac) //user supplied gradient
                calcG(&ndec,&ndat,x,fjac,&grad);     
            else {//internal djacobix
                if(djacobix(&calcNDG,&ndec,&ndat,fjac,x,&tol,&fun) != TR_SUCCESS)
                {
                    *exitflag = -6;
                    MKL_Free_Buffers();
                    return;
                }
            }
                
        }
    }
    
    //Get Solution Information
    if(con)
        res = dtrnlspbc_get (&handle, &miter, &st_cr, &r1, &r2);
    else
        res = dtrnlsp_get (&handle, &miter, &st_cr, &r1, &r2);
        
    if (res != TR_SUCCESS)
	{
        *exitflag = -7;
        MKL_Free_Buffers();
        return;
	}
    
    //Calculate final residual
    *fval = 0;
    for(i = 0; i < ndat; i++)
    	*fval += fvec[i]*fvec[i];   
    
    //Check if maxtime exceeded
    if(((double)(end-start))/CLOCKS_PER_SEC > maxtime)
        st_cr = 15;
    
    //Assign solution
    *exitflag = getStatus(st_cr);
    *iter = miter;
    *feval = (double)uiparm[0];
    *geval = (double)jacEval;

	// Free handle memory
    if(con)
        res = dtrnlspbc_delete (&handle);
    else
        res = dtrnlsp_delete (&handle);
        
	if (res != TR_SUCCESS)
	{
		*exitflag = -8;
        MKL_Free_Buffers();
        return;
	}

    //Free allocated memory
	free (fvec);
	free (fjac);                 
	MKL_Free_Buffers();
    if(lb) mxFree(lb);
    if(ub) mxFree(ub);

    //Print Header
    if(uiparm[1]){            
        //Termination Detected
        switch(st_cr)
        {
            //Success
            case 2:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** trust region convergence: delta < %g ***\n",eps[0]); break;
            case 3:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** function convergence: ||F(x)||2 < %g ***\n",eps[1]); break;
            case 5:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** trial step convergence: ||s||2 < %g ***\n",eps[3]); break;
            case 6:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** ||F(x)||2 - ||F(x) - J(x)s||2 < %g ***\n",eps[4]); break;
            //Error
            case 1:
                mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;  
            case 15:
                mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;  
            //Early Exit
            case 4:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: the Jacobian matrix is singular. ***\n"); break;
            //Other Error
            default:
                mexPrintf("\n *** ERROR: internal error code %d ***\n",(int)st_cr); break;
        }
        
        if(*exitflag==1)
            mexPrintf("\n Final SSE: %12.5g\n In %3.0f iterations\n",*fval,*iter);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
}

//Function Value Calculation
int calcF(MKL_INT *ndec, MKL_INT *ndat, double *x, double *f, double *ydata, void *userdata, int *uiparm)
{
    bool stop = false, havrnorm = false;
    MKL_INT i, stat;
    double *fval, rnorm, evaltime;   

    user_function_data *fun = (user_function_data*)userdata;
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting MKLTRNLS...\n\n");
        return -1; //terminate
    }
    
    //Check for maxtime expiry    
    if(evaltime > maxtime)
    {
        mexPrintf("\nMaximum Solver Time Exceeded. Exiting MKLTRNLS...\n\n");
        return -1; //terminate
    }
    
    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, *ndec * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");
    
    //Get Objective
    fval = mxGetPr(fun->plhs[0]);

    //Assign Objective + subtract ydata
    for(i=0;i<*ndat;i++)
        f[i] = fval[i]-ydata[i];
    
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
    
    //Iteration Printing
    if(uiparm[1] > 1) {       
        //Calculate Residual Norm
        rnorm = 0; havrnorm = true;
        for(i=0;i<*ndat;i++)
            rnorm += f[i]*f[i];
        
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
            for(i=0;i<*ndat;i++)
                rnorm += f[i]*f[i];

        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), uiparm , sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), &rnorm, sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), x, *ndec * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");
        
        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
    }    
    //Increment feval
    uiparm[0]++;
    
    //Check for iterfun terminate
    if (stop) {
        mexPrintf("\nIterFun called Stop. Exiting MKLTRNLS...\n\n");
        return -1; //terminate
    }
    
    return 1;
}

//Function Gradient Calculation
void calcG(MKL_INT *ndec, MKL_INT *ndat, double *x, double *g, void *userdata)
{
    int stat;
    double *grad;   
    
    user_function_data *fun = (user_function_data*)userdata;

    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, *ndec * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
      mexErrMsgTxt("Error calling Gradient Function!");
    
    //Get Gradient
    grad = mxGetPr(fun->plhs[0]);
    //Assign Gradient
    memcpy(g,grad,*ndec**ndat*sizeof(double));
    
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
}

//DJACOBI Callback Calculation
void calcNDG(MKL_INT *m, MKL_INT *n, double *x, double *f, void *data)
{
    int stat;

    user_function_data *fun = (user_function_data*)data;
    
    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, *n * sizeof(double));
    
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");
    
    //Check Return Size
    if(mxGetNumberOfElements(fun->plhs[0]) != *m)
        mexErrMsgTxt("Incorrect sized vector returned from function");
    
    //Save Result
    memcpy(f, mxGetPr(fun->plhs[0]), *m * sizeof(double));

    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
}

void checkInputs(const mxArray *prhs[], int nrhs, int *con)
{    
    size_t Mx0;
    
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 4 arguments to mklTRnls!\n\nmklTRnls(fun,grad,x0,ydata)\n");
       
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
            *con = 1;
    }
    
    //Check Options
    if(nrhs > 6) {
        if(!mxIsStruct(prhs[6]))
            mexErrMsgTxt("The specified options must be a structure!");
    }

}

double getStatus(int stat)
{
    switch((int)stat)
    {     
        case 2:         //stopped by TR
        case 3:         //stopped by ftol
        case 5:         //stopped by step size
        case 6:         //stopped by ftol - f-gs
            return 1;
            break;
        case 1:         //itmax
        case 15:        //max time
            return 0;
            break;
        case 4:         //singular
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
    mexPrintf(" MKLTRNLS: Intel MKL Trust Region Nonlinear Least Squares [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);  
    PRINT_BUILD_INFO;
    mexPrintf("  - Released as part of the Intel Math Kernel Library\n  - http://software.intel.com/en-us/intel-mkl\n");
    
    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}