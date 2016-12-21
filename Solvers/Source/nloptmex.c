/* NLOPTMEX - A MATLAB MEX Interface to NLOPT
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.inverseproblem.co.nz
 */

/* Based in parts on nlopt_optimize-mex.c supplied with NLOPT */

#include <mex.h>
#include <math.h>
#include <time.h>

#include "nlopt.h"
#include "config.h"
#include "opti_util.h"

#define FLEN 256 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */    
//Iteration callback structure
typedef struct {
    char f[FLEN];
    mxArray *plhs[1];
    mxArray *prhs[3];
    bool enabled;
} iter_fun_data;    
//User Structures
typedef struct {
     char f[FLEN], g[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     mxArray *prhs_g[MAXRHS];
     int xrhs, nrhs, xrhs_g, nrhs_g;
     int nfeval, ngeval, printLevel;
     clock_t start;
     iter_fun_data *iterF;
     nlopt_opt *opt;
} obj_fun_data;
typedef struct {
     char f[FLEN], g[FLEN];
     mxArray *plhs[1];
     mxArray *prhs[MAXRHS];
     mxArray *prhs_g[MAXRHS];
     int xrhs, nrhs, xrhs_g, nrhs_g;
     char *contype, mode;
     double *nlrhs, *nle;
     int ncval, njval, ncon;
     nlopt_opt *opt;
} con_fun_data;

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
bool checkIdentX(const double *x, double *X, int n);
static double obj_callback(unsigned int n, const double *x, double *gradient, void *user_data);
static void con_callback(unsigned m, double *result, unsigned n, const double *x, double *gradient, void *user_data);

//Macros
#define CHECK(cond, msg) if (!(cond)) { mexErrMsgTxt(msg); }
#define NLOPT_ERR(cond, msg) if (!(cond)) { nlopt_destroy(opt); nlopt_destroy(local_opt); mexErrMsgTxt(msg); };

//Argument Enums (in expected order of arguments)
enum {eNLPROB, eX0};                   
//PRHS Defines    
#define pNLPROB     prhs[eNLPROB]
#define pX0         prhs[eX0]
//Field Defines
#define pFUN        mxGetField(pNLPROB,0,"objective")
#define pGRAD       mxGetField(pNLPROB,0,"gradient")
#define pNLCON      mxGetField(pNLPROB,0,"nlcon")
#define pNLJAC      mxGetField(pNLPROB,0,"nljac")
#define pNLRHS      mxGetField(pNLPROB,0,"nlrhs")
#define pNLE        mxGetField(pNLPROB,0,"nle")
#define pLB         mxGetField(pNLPROB,0,"lb")
#define pUB         mxGetField(pNLPROB,0,"ub")
#define pOPTS       mxGetField(pNLPROB,0,"options")
#define pLOPTS      mxGetField(mxGetField(pNLPROB,0,"options"),0,"local_optimizer")

#define INEQ -1
#define EQ 1

//Ctrl-C Detection
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

//REQUIRED FOR LINKER [seems odd..]
nlopt_algorithm nlopt_local_search_alg_deriv;
nlopt_algorithm nlopt_local_search_alg_nonderiv;
int nlopt_local_search_maxeval;
unsigned nlopt_stochastic_population;

//Global X, C, J for avoiding multiple constraint and jacobian evaluations
double *X = NULL, *C = NULL, *J = NULL;

//Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *lb = NULL, *ub = NULL, *x0 = NULL;
    
    //Return Args
    double *x, *fval, *exitflag, *fvals, *gvals, *cvals, *jvals;
    const char *fnames[4] = {"objective","constraint","gradient","jacobian"};
    
    //Options    
    int maxfeval = 1500, lmaxfeval = 1500;
    int printLevel = 0;
    unsigned pop = 0, lpop = 0, vecstore = 0;
    double maxtime = 1000, lmaxtime = 1000;
    double tolrfun = 1e-7, ltolrfun = 1e-7;
    double tolafun = 1e-7, ltolafun = 1e-7;    
    double stopval = -HUGE_VAL;
    
    //NLOPT Vars
    nlopt_opt opt = NULL, local_opt = NULL;
    nlopt_algorithm algorithm = NLOPT_NUM_ALGORITHMS;
    nlopt_algorithm lalgorithm = NLOPT_NUM_ALGORITHMS;
    nlopt_result ret;
    
    //Internal Vars
    size_t ndec = 0, nineq = 0, neq = 0, ncon = 0, i;
    double *alb, *aub, *zeros, *intol = NULL, *eqtol = NULL, *nle;
    char *contype = NULL;
    
    //User Function Structure
    obj_fun_data fun;
    con_fun_data icon;
    con_fun_data econ;
    //Iteration Callback Structure
    iter_fun_data iterF;
    
    //Check for ver display
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(PACKAGE_VERSION);
        return;
    } 
    
    //Check Inputs
    checkInputs(prhs,nrhs); 
    
    //Get Sizes
    ndec = mxGetNumberOfElements(pX0);
    //Get Objective Function Handle
    if (mxIsChar(pFUN)) {
        CHECK(mxGetString(pFUN, fun.f, FLEN) == 0,"Error reading Objective function");
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
    if(pGRAD && !mxIsEmpty(pGRAD)) {
        if (mxIsChar(pGRAD)) {
            CHECK(mxGetString(pGRAD, fun.g, FLEN) == 0,"Error reading Gradient function");
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
    else
        fun.nrhs_g = 0;
    //Get Constraints Handles (if exists)
    if(pNLCON && !mxIsEmpty(pNLCON) && !mxIsEmpty(pNLE)) {
        //Determine number of ineq and eq
        ncon = mxGetNumberOfElements(pNLE);
        nle = mxGetPr(pNLE);
        for(i = 0; i < ncon; i++)
            if(nle[i] == 0)
                neq++;
            else
                nineq++;     
        //Based on above, save to respective structure
        if (mxIsChar(pNLCON)) {
            if(nineq) {
                CHECK(mxGetString(pNLCON, icon.f, FLEN) == 0,"error reading constraint name string");
                icon.nrhs = 1; icon.xrhs = 0;
            }
            if(neq) {
                CHECK(mxGetString(pNLCON, econ.f, FLEN) == 0,"error reading constraint name string");
                econ.nrhs = 1; econ.xrhs = 0;
            }
        } else {
            if(nineq) {
                icon.prhs[0] = (mxArray*)pNLCON;
                strcpy(icon.f, "feval");
                icon.nrhs = 2; icon.xrhs = 1;
            }
            if(neq) {
                econ.prhs[0] = (mxArray*)pNLCON;
                strcpy(econ.f, "feval");
                econ.nrhs = 2; econ.xrhs = 1;
            }
        }
        //Create x and set nlrhs, nle for each
        if(nineq) {
            icon.prhs[icon.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL);
            icon.nlrhs = mxGetPr(pNLRHS);
            icon.nle = mxGetPr(pNLE);
            icon.ncon = (int)mxGetNumberOfElements(pNLE);                            
        }
        if(neq) {
            econ.prhs[econ.xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL);
            econ.nlrhs = mxGetPr(pNLRHS);
            econ.nle = mxGetPr(pNLE);
            econ.ncon = (int)mxGetNumberOfElements(pNLE);
        }
        //Create contype array
        contype = mxCalloc(ncon,sizeof(char));
        for(i=0;i<ncon;i++)
            if(nle[i] == 0)
                contype[i] = EQ;
            else
                contype[i] = INEQ;
        //Assign to structures
        icon.contype = contype;
        econ.contype = contype;
    }
    //Get Jacobian Function Handle (if exists)
    if(pNLJAC && !mxIsEmpty(pNLJAC)) {
        //Based on above, save to respective structure        
        if (mxIsChar(pNLJAC)) {
            if(nineq) {
                CHECK(mxGetString(pNLJAC, icon.g, FLEN) == 0,"error reading jacobian name string");
                icon.nrhs_g = 1; icon.xrhs_g = 0;
            }
            if(neq) {
                CHECK(mxGetString(pNLJAC, econ.g, FLEN) == 0,"error reading jacobian name string");
                econ.nrhs_g = 1; econ.xrhs_g = 0;
            }
        } else {
            if(nineq) {
                icon.prhs_g[0] = (mxArray*)pNLJAC;
                strcpy(icon.g, "feval");
                icon.nrhs_g = 2; icon.xrhs_g = 1;
            }
            if(neq) {
                econ.prhs_g[0] = (mxArray*)pNLJAC;
                strcpy(econ.g, "feval");
                econ.nrhs_g = 2; econ.xrhs_g = 1;
            }
        }   
        //Create x for each
        if(nineq) icon.prhs_g[icon.xrhs_g] = mxCreateDoubleMatrix(ndec, 1, mxREAL);   
        if(neq)   econ.prhs_g[econ.xrhs_g] = mxCreateDoubleMatrix(ndec, 1, mxREAL);
    }
    else {
        icon.nrhs_g = 0;
        econ.nrhs_g = 0;
    }
    
    //Get x0
    x0 = mxGetPr(pX0);
    
    //Get Options if specified
    if(pOPTS && !mxIsEmpty(pOPTS)) {
        if(mxGetField(pOPTS,0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display"));
        if(mxGetField(pOPTS,0,"maxfeval"))
            maxfeval = (int)*mxGetPr(mxGetField(pOPTS,0,"maxfeval"));
        if(mxGetField(pOPTS,0,"maxtime"))
            maxtime = *mxGetPr(mxGetField(pOPTS,0,"maxtime"));
        if(mxGetField(pOPTS,0,"tolafun"))
            tolafun = *mxGetPr(mxGetField(pOPTS,0,"tolafun"));
        if(mxGetField(pOPTS,0,"tolrfun"))
            tolrfun = *mxGetPr(mxGetField(pOPTS,0,"tolrfun")); 
        if(mxGetField(pOPTS,0,"iterfun") && !mxIsEmpty(mxGetField(pOPTS,0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(pOPTS,0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
            //Set structure
            fun.iterF = &iterF;
        }
        else
            fun.iterF = NULL;
        
        //Optimizer Config
        if(mxGetField(pOPTS,0,"algorithm"))
            algorithm = (nlopt_algorithm)*mxGetPr(mxGetField(pOPTS,0,"algorithm"));
        if(mxGetField(pOPTS,0,"initial_pop") && !mxIsEmpty(mxGetField(pOPTS,0,"initial_pop")))
            pop = (int)*mxGetPr(mxGetField(pOPTS,0,"initial_pop"));
        if(mxGetField(pOPTS,0,"vector_storage") && !mxIsEmpty(mxGetField(pOPTS,0,"vector_storage")))
            vecstore = (int)*mxGetPr(mxGetField(pOPTS,0,"vector_storage"));
        if(mxGetField(pOPTS,0,"stopval") && !mxIsEmpty(mxGetField(pOPTS,0,"stopval")) && !mxIsInf(*mxGetPr(mxGetField(pOPTS,0,"stopval"))))
            stopval = *mxGetPr(mxGetField(pOPTS,0,"stopval"));      
    }   

    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateStructMatrix(1,1,4,fnames);
    mxSetField(plhs[3],0,fnames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[2],mxCreateDoubleMatrix(1,1, mxREAL)); 
    mxSetField(plhs[3],0,fnames[3],mxCreateDoubleMatrix(1,1, mxREAL)); 
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);     
    fvals  = mxGetPr(mxGetField(plhs[3],0,fnames[0]));
    cvals  = mxGetPr(mxGetField(plhs[3],0,fnames[1]));
    gvals  = mxGetPr(mxGetField(plhs[3],0,fnames[2]));
    jvals  = mxGetPr(mxGetField(plhs[3],0,fnames[3]));    
    
    //Make NLOPT Object  
    NLOPT_ERR( opt = nlopt_create(algorithm, (int)ndec), "Error creating NLOPT structure. Ensure algorithm is set correctly." );
    
    //Set Bounds
    alb = mxCalloc(ndec,sizeof(double));
    aub = mxCalloc(ndec,sizeof(double));
    if(pLB && !mxIsEmpty(pLB)) {
        lb = mxGetPr(pLB);
        for(i = 0; i < (int)ndec; i++)
            if(mxIsInf(lb[i]))
                alb[i] = -HUGE_VAL;
            else
                alb[i] = lb[i];
    }        
    else
        for(i = 0; i < ndec; i++)
            alb[i] = -HUGE_VAL;
    
    if(pUB && !mxIsEmpty(pUB)) {
        ub = mxGetPr(pUB);
        for(i = 0; i < (int)ndec; i++)
            if(mxIsInf(ub[i]))
                aub[i] = HUGE_VAL;
            else
                aub[i] = ub[i];
    }        
    else 
        for(i = 0; i < ndec; i++)
            aub[i] = HUGE_VAL;
    
    NLOPT_ERR( nlopt_set_lower_bounds(opt, alb), "Error setting lower bounds" );
    NLOPT_ERR( nlopt_set_upper_bounds(opt, aub), "Error setting upper bounds" );
    
    //Set User Options
    NLOPT_ERR( nlopt_set_ftol_rel(opt, tolrfun), "Error setting ftol_rel" );
    NLOPT_ERR( nlopt_set_ftol_abs(opt, tolafun), "Error setting ftol_abs" );
    NLOPT_ERR( nlopt_set_maxeval(opt, maxfeval), "Error setting maxeval" );
    NLOPT_ERR( nlopt_set_maxtime(opt, maxtime), "Error setting maxtime" );
    NLOPT_ERR( nlopt_set_stopval(opt, stopval), "Error setting stopval" );
    NLOPT_ERR( nlopt_set_population(opt, pop), "Error setting population" );
    NLOPT_ERR( nlopt_set_vector_storage(opt, vecstore), "Error setting vector storage" );
    
    //Set Default Options
    zeros = mxCalloc(ndec,sizeof(double));
    for(i = 0; i < (int)ndec; i++)
        zeros[i] = 0.0;
    
    NLOPT_ERR( nlopt_set_xtol_rel(opt, 0), "Error setting xtol_rel" );
    NLOPT_ERR( nlopt_set_xtol_abs(opt, zeros ), "Error setting xtol_abs" );        
    
    //Check for local optimizer options
    if(pOPTS && !mxIsEmpty(pOPTS) && pLOPTS && !mxIsEmpty(pLOPTS)) {
        if(mxGetField(pLOPTS,0,"algorithm"))
            lalgorithm = (nlopt_algorithm)*mxGetPr(mxGetField(pLOPTS,0,"algorithm"));
        //Set Algorithm
        NLOPT_ERR( local_opt = nlopt_create(lalgorithm, (int)ndec), "Error creating NLOPT local optimizer structure. Ensure local_optimizer.algorithm is set correctly." );
        //Get local options
        if(mxGetField(pLOPTS,0,"maxfeval"))
            lmaxfeval = (int)*mxGetPr(mxGetField(pLOPTS,0,"maxfeval"));
        if(mxGetField(pLOPTS,0,"maxtime"))
            lmaxtime = *mxGetPr(mxGetField(pLOPTS,0,"maxtime"));
        if(mxGetField(pLOPTS,0,"tolrfun"))
            ltolrfun = *mxGetPr(mxGetField(pLOPTS,0,"tolrfun"));
        if(mxGetField(pLOPTS,0,"tolafun"))
            ltolafun = *mxGetPr(mxGetField(pLOPTS,0,"tolafun"));
        if(mxGetField(pLOPTS,0,"initial_pop"))
            lpop = (int)*mxGetPr(mxGetField(pLOPTS,0,"initial_pop"));
        //Set local options
        NLOPT_ERR( nlopt_set_ftol_rel(local_opt, ltolrfun), "Error setting local ftol_rel" );
        NLOPT_ERR( nlopt_set_ftol_abs(local_opt, ltolafun), "Error setting local ftol_abs" );
        NLOPT_ERR( nlopt_set_maxeval(local_opt, lmaxfeval), "Error setting local maxfeval" );
        NLOPT_ERR( nlopt_set_maxtime(local_opt, lmaxtime), "Error setting local maxtime" );   
        NLOPT_ERR( nlopt_set_population(local_opt, lpop), "Error setting local population" );
        //Set default options
        NLOPT_ERR( nlopt_set_stopval(local_opt, -HUGE_VAL), "Error setting local stop val");
        NLOPT_ERR( nlopt_set_xtol_rel(local_opt, 0), "Error setting local xtol_rel" );
        NLOPT_ERR( nlopt_set_xtol_abs(local_opt, zeros ), "Error setting local xtol_abs" );        
        NLOPT_ERR( nlopt_set_vector_storage(local_opt, 0), "Error setting local vector storage" );
        //Set Local Optimizer
        NLOPT_ERR( nlopt_set_local_optimizer(opt, local_opt), "Error setting NLOPT local optimizer" );
        nlopt_destroy(local_opt); local_opt = NULL;
    }
    
    //Set Objective & Constraints
    NLOPT_ERR( nlopt_set_min_objective(opt, obj_callback, &fun), "Error setting min objective");
    if(nineq > 0) {
        intol = (double*)mxCalloc(nineq,sizeof(double));
        for(i = 0; i < (int)nineq; i++)
            intol[i] = 1e-8;
        NLOPT_ERR( nlopt_add_inequality_mconstraint(opt, (int)nineq, con_callback, &icon, intol),  "Error adding Inequality constraints");
    }
    if(neq > 0) {
        eqtol = (double*)mxCalloc(neq,sizeof(double));
        for(i = 0; i < (int)neq; i++)
            eqtol[i] = 1e-8;
        NLOPT_ERR( nlopt_add_equality_mconstraint(opt, (int)neq, con_callback, &econ, eqtol),  "Error adding Equality constraints");
    }
    
    //Set Print Level
    fun.printLevel = printLevel;
    
    //Set Default Counters
    fun.nfeval = fun.ngeval = 0;
    icon.ncval = icon.njval = 0;
    econ.ncval = econ.njval = 0;
    
    //Copy in NLOPT Object to User Structures
    fun.opt = &opt;
    icon.opt = &opt; icon.mode = INEQ;
    econ.opt = &opt; econ.mode = EQ;
    
    //Copy initial guess to x
    memcpy(x,x0,ndec*sizeof(double)); 
    
    //Ensure x within problem bounds
    for(i=0;i<(int)ndec;i++)
    {
        if(x[i] < alb[i])
            x[i] = alb[i];
        else if(x[i] > aub[i])
            x[i] = aub[i];
    }
    
    //Create Global Memory for storing previous evaluations
    X = mxCalloc(ndec,sizeof(double));
    if(ncon) {
        C = mxCalloc(ncon,sizeof(double));
        if(icon.nrhs_g || econ.nrhs_g)
            J = mxCalloc(ncon*ndec,sizeof(double));
    }
    
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is NLOPT v%s\n",PACKAGE_VERSION);  
        mexPrintf(" Author: Steven G. Johnson\n\n");
        mexPrintf(" Algorithm: %s\n",nlopt_algorithm_name(algorithm));
        if(pLOPTS && lalgorithm != NLOPT_NUM_ALGORITHMS)
            mexPrintf(" Local Algorithm: %s\n\n",nlopt_algorithm_name(lalgorithm));
        else
            mexPrintf("\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Inequality Constraints: %4d\n",nineq);
        mexPrintf(" # Equality Constraints:   %4d\n",neq);

        mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;");
    }
 
    //Start Timer    
    fun.start = clock();

    //Solve Problem
    ret = nlopt_optimize(opt, x, fval);
        
    //Set Outputs
    *exitflag = (double)ret;
    *fvals = (double)fun.nfeval;
    *gvals = (double)fun.ngeval;
    *cvals = (double)icon.ncval + (double)econ.ncval;
    *jvals = (double)icon.njval + (double)econ.njval;
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch((int)*exitflag)
        {
            //Success
            case 1:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Generic Success ***\n"); break;
            case 2:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Optimization stopped because stopval was reached ***\n"); break;
            case 3:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Optimization stopped because ftol_rel or ftol_abs was reached ***\n"); break;
            case 4:
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n *** Optimization stopped because xtol_rol or xtol_abs was reached ***\n"); break;
            //Error
            case 5:
                mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); break;                      
            case 6:
                mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;  
            case -1:
                mexPrintf("\n *** ERROR: Generic Failure ***\n"); break;
            case -2:
                mexPrintf("\n *** ERROR: Invalid Arguments ***\n"); break;    
            case -3:
                mexPrintf("\n *** ERROR: Ran out of Memory ***\n"); break;       
            case -4:
                mexPrintf("\n *** ERROR: Halted because of roundoff errors ***\n"); break;
            //Early Exit
            case -5:
                mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
            //Other Error
            default:
                mexPrintf("\n *** ERROR: internal error code %d ***\n",(int)*exitflag); break;
        }
        
        if(*exitflag==1)
            mexPrintf("\n Final fval: %12.5g\n In %3.0f function evaluations\n",*fval,*fvals);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Clean Up
    nlopt_destroy(opt);
    mxFree(alb); alb = NULL;
    mxFree(aub); aub = NULL;
    mxFree(zeros); zeros = NULL;
    if(intol) {mxFree(intol); intol = NULL;}
    if(eqtol) {mxFree(eqtol); eqtol = NULL;}
    if(contype) {mxFree(contype); contype = NULL;}
    if(X) {mxFree(X); X = NULL;}
    if(C) {mxFree(C); C = NULL;}
    if(J) {mxFree(J); J = NULL;}
}

static double obj_callback(unsigned int n, const double *x, double *gradient, void *user_data)
{
    bool stop = false;
    int stat;
    obj_fun_data *fun = (obj_fun_data *) user_data;
    double f, *g, evaltime;
    clock_t end;

    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-fun->start))/CLOCKS_PER_SEC;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting NLOPT...\n\n");
        nlopt_force_stop(*fun->opt);
        return 0;
    }  
    
    //Set x
    fun->plhs[0] = NULL;
    memcpy(mxGetPr(fun->prhs[fun->xrhs]), x, n * sizeof(double));
    //Call MATLAB
    stat = mexCallMATLAB(1, fun->plhs, fun->nrhs, fun->prhs, fun->f);
    if(stat)
        mexErrMsgTxt("Error calling Objective Function!");
    //Get Objective
    f = *mxGetPr(fun->plhs[0]);
    // Clean up Ptr
    mxDestroyArray(fun->plhs[0]);
    //Increment counter
    fun->nfeval++;

    //Optionally Get Gradient
    if(gradient) {    
        //Check we have a gradient function
        if(fun->nrhs_g == 0)
            mexErrMsgTxt("The selected algorithm requires a gradient, but you have not supplied one!");
        //Get Gradient
        fun->plhs[0] = NULL;
        memcpy(mxGetPr(fun->prhs_g[fun->xrhs_g]), x, n * sizeof(double));

        stat = mexCallMATLAB(1, fun->plhs, fun->nrhs_g, fun->prhs_g, fun->g);
        if(stat)
            mexErrMsgTxt("Error calling Gradient Function!");
        if(mxGetNumberOfElements(fun->plhs[0]) != n)
            mexErrMsgTxt("Incorrect number of elements returned from Gradient");
        if(mxIsSparse(fun->plhs[0]) || !mxIsNumeric(fun->plhs[0]) || mxIsComplex(fun->plhs[0]))
            mexErrMsgTxt("The Gradient must return a real dense vector");
        //Get Gradient
        g = mxGetPr(fun->plhs[0]);
        //Assign Gradient
        memcpy(gradient,g,n*sizeof(double));
        // Clean up Ptr
        mxDestroyArray(fun->plhs[0]);
        //Increment counter
        fun->ngeval++;
    }

    //Iteration Printing
    if(fun->printLevel > 1) {               
        if(fun->nfeval == 1 || !(fun->nfeval%10))
            mexPrintf(" feval       time           fval\n");

        mexPrintf("%5d       %5.2f    %12.5g\n",fun->nfeval,evaltime,f);
        mexEvalString("drawnow;"); //flush draw buffer
    }

    //Iteration Callback
    if(fun->iterF != NULL)
    {
        fun->iterF->plhs[0] = NULL;
        memcpy(mxGetData(fun->iterF->prhs[1]), &fun->nfeval, sizeof(int));
        memcpy(mxGetPr(fun->iterF->prhs[2]), &f, sizeof(double));
        memcpy(mxGetPr(fun->iterF->prhs[3]), x, n * sizeof(double));
        stat = mexCallMATLAB(1, fun->iterF->plhs, 4, fun->iterF->prhs, fun->iterF->f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(fun->iterF->plhs[0]);
        // Clean up Ptr
        mxDestroyArray(fun->iterF->plhs[0]);

        if(stop) {
            mexPrintf("\nIterFun called Stop. Exiting NLOPT...\n\n");
            nlopt_force_stop(*fun->opt);
        }
    }
    return f;
}
    
static void con_callback(unsigned m, double *result, unsigned n, const double *x, double *gradient, void *user_data)
{
    unsigned int i, j, k=0, stat;
    con_fun_data *con = (con_fun_data*)user_data;
    double *c, *g;

    //Check for identical X (in which case we can avoid re-evaluating constraint and jacobian functions between ineq and eq calls)
    if(con->ncval<1 || !checkIdentX(x,X,n)) {  
        //Set x
        con->plhs[0] = NULL;
        memcpy(mxGetPr(con->prhs[con->xrhs]), x, n * sizeof(double));
        //Copy to saved X
        memcpy(X,x,n*sizeof(double));
        //Call MATLAB
        stat = mexCallMATLAB(1, con->plhs, con->nrhs, con->prhs, con->f);    
        if(stat)
            mexErrMsgTxt("Error calling Constraint Function!");
        //Check result
        if(mxGetNumberOfElements(con->plhs[0]) != con->ncon)
            mexErrMsgTxt("Incorrect number of elements returned from constraint function");
        if(mxIsSparse(con->plhs[0]) || !mxIsNumeric(con->plhs[0]) || mxIsComplex(con->plhs[0]))
            mexErrMsgTxt("Constraint function must return a real dense vector");

        //Get Result
        c = mxGetPr(con->plhs[0]);
        //Save C
        memcpy(C,c,con->ncon*sizeof(double));
        //Index Elements of this constraint (ineq or eq)
        for(i=0;i<(size_t)con->ncon;i++) {
            if(con->contype[i] == con->mode) {
                if(con->nle[i] == 1)
                    result[k++] = -c[i] + con->nlrhs[i];
                else
                    result[k++] = c[i] - con->nlrhs[i];
            }   
        }    
        // Clean up Ptr
        mxDestroyArray(con->plhs[0]);
        //Increment counter
        con->ncval++;        
        
        //Optionally Get Jacobian
        if(gradient) {    
            //Check we have a jacobian function
            if(con->nrhs_g == 0)
                mexErrMsgTxt("The selected algorithm requires a Jacobian, but you have not supplied one!");
            //Get Jacobian
            con->plhs[0] = NULL;
            memcpy(mxGetPr(con->prhs_g[con->xrhs_g]), x, n * sizeof(double));
            //Call MATLAB
            stat = mexCallMATLAB(1, con->plhs, con->nrhs_g, con->prhs_g, con->g);
            if(stat)
                mexErrMsgTxt("Error calling Jacobian Function!");
            if(mxGetNumberOfElements(con->plhs[0]) != con->ncon*n)
                mexErrMsgTxt("Incorrect number of elements returned from Jacobian");
            if(mxIsSparse(con->plhs[0]) || !mxIsNumeric(con->plhs[0]) || mxIsComplex(con->plhs[0]))
                mexErrMsgTxt("The Jacobian must return a real dense vector");
            //Get Jacobian
            g = mxGetPr(con->plhs[0]);
            //Save J
            memcpy(J,g,con->ncon*n*sizeof(double));
            //Assign Jaobian, transposing as we go
            k = 0;
            for(i=0;i<(size_t)con->ncon;i++)
                if(con->contype[i] == con->mode)
                    for(j=0;j<n;j++)
                        gradient[k++] = g[i + j*con->ncon];            
            // Clean up Ptr
            mxDestroyArray(con->plhs[0]);
            //Increment counter
            con->njval++;            
        }
    }
    //Else is indentical
    else {
        k = 0;
        //Process Saved Constraints to result arg
        for(i=0;i<(size_t)con->ncon;i++) {
            if(con->contype[i] == con->mode) {
                if(con->nle[i] == 1)
                    result[k++] = -C[i] + con->nlrhs[i];
                else
                    result[k++] = C[i] - con->nlrhs[i];
            }   
        } 
        //Optionally process saved jacobian to gradient arg
        if(gradient){
            k = 0;
            for(i=0;i<(size_t)con->ncon;i++)
                if(con->contype[i] == con->mode)
                    for(j=0;j<n;j++)
                        gradient[k++] = J[i + j*con->ncon]; 
        }
    }    
}


void checkInputs(const mxArray *prhs[], int nrhs)
{   
    size_t ndec, i;
    double *nle;
    
    if(nrhs < 2)
        mexErrMsgTxt("You must supply at least 2 arguments to nlopt!\n\nlopt(nlprob,x0)\n");
       
    //Check Types
    if(!mxIsStruct(pNLPROB) || mxIsEmpty(pNLPROB))
        mexErrMsgTxt("nlprob must be a structure!");
    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");
    if(pOPTS && !mxIsStruct(pOPTS))
    	mexErrMsgTxt("The specified options must be a structure!");

    //Check we have an algorithm
    if(pOPTS == NULL || mxIsEmpty(pOPTS) || mxGetField(pOPTS,0,"algorithm") == NULL || mxIsEmpty(mxGetField(pOPTS,0,"algorithm")))
        mexErrMsgTxt("You must supply a field nlprob.options.algorithm which specifies the NLOPT algorithm number");
    if(*mxGetPr(mxGetField(pOPTS,0,"algorithm")) < 0 || *mxGetPr(mxGetField(pOPTS,0,"algorithm")) >= NLOPT_NUM_ALGORITHMS)
        mexErrMsgTxt("NLOPT algorithm number is not valid");
    //Check local algorithm if specified
    if(pLOPTS != NULL && !mxIsEmpty(pLOPTS)) {
        if(!mxIsStruct(pLOPTS))
            mexErrMsgTxt("nlprob.options.local_optimizer must be a structure!");
        if(mxGetField(pLOPTS,0,"algorithm") == NULL || mxIsEmpty(mxGetField(pLOPTS,0,"algorithm")))
            mexErrMsgTxt("You must supply a field nlprob.options.local_optimizer.algorithm which specifies the NLOPT algorithm number when using a local optimizer");
        if(*mxGetPr(mxGetField(pLOPTS,0,"algorithm")) < 0 || *mxGetPr(mxGetField(pLOPTS,0,"algorithm")) >= NLOPT_NUM_ALGORITHMS)
            mexErrMsgTxt("NLOPT local optimizer algorithm number is not valid");
    }
        
    //Check inputs
    ndec = mxGetNumberOfElements(pX0);
    //Objective & Gradient
    if(pFUN == NULL || mxIsEmpty(pFUN) || (!mxIsChar(pFUN) && !mxIsFunctionHandle(pFUN)))
        mexErrMsgTxt("nlprob.objective must be a function name (string) or function handle!");
    if(pGRAD && !mxIsEmpty(pGRAD) && !mxIsChar(pGRAD) && !mxIsFunctionHandle(pGRAD))
        mexErrMsgTxt("nlprob.gradient must be a function name (string) or function handle!");
    //Bounds
    if(pLB && !mxIsEmpty(pLB) && (mxGetNumberOfElements(pLB) != ndec || mxIsComplex(pLB)))
        mexErrMsgTxt("nlprob.lb must be a real column vector with the same number of elements as x0");
    if(pUB && !mxIsEmpty(pUB) && (mxGetNumberOfElements(pUB) != ndec || mxIsComplex(pUB)))
        mexErrMsgTxt("nlprob.ub must be a real column vector with the same number of elements as x0");
    //Constraints
    if(pNLCON && !mxIsEmpty(pNLCON)) {
        if(!mxIsChar(pNLCON) && !mxIsFunctionHandle(pNLCON))
            mexErrMsgTxt("nlprob.nlcon must be a function name (string) or function handle!");
        if(pNLRHS == NULL || mxIsEmpty(pNLRHS))
            mexErrMsgTxt("You must supply nlprob.nlrhs when specifying nonlinear constraints");
        if(pNLE == NULL || mxIsEmpty(pNLE))
            mexErrMsgTxt("You must supply nlprob.nle when specifying nonlinear constraints");          
        if(mxGetNumberOfElements(pNLRHS) != mxGetNumberOfElements(pNLE))
            mexErrMsgTxt("nlprob.nlrhs and nlprob.nle are not the same length!");
        if(mxIsSparse(pNLRHS) || mxIsComplex(pNLRHS) || !mxIsDouble(pNLRHS))
            mexErrMsgTxt("nlprob.nlrhs must be a real dense column vector");
        if(mxIsSparse(pNLE) || mxIsComplex(pNLE) || !mxIsDouble(pNLE))
            mexErrMsgTxt("nlprob.nle must be a real dense column vector");
        nle = mxGetPr(pNLE);
        for(i=0;i<mxGetNumberOfElements(pNLE);i++)
            if(nle[i] < -1 || nle[i] > 1)
                mexErrMsgTxt("nlprob.nle must only contain -1, 0 or 1!");
    }   
    //Jacobian
    if(pNLJAC && !mxIsEmpty(pNLJAC) && !mxIsChar(pNLJAC) && !mxIsFunctionHandle(pNLJAC))
        mexErrMsgTxt("nlprob.nljac must be a function name (string) or function handle!");
}

bool checkIdentX(const double *x, double *X, int n)
{
    int i;
    //Check each for equality, should be no need for tolerance
    for(i=0;i<n;i++)
        if(x[i]!=X[i])
            return false;
    return true;
}

//Print Solver Information
void printSolverInfo()
{    
    char vbuf[6]; getVSVer(vbuf);    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" NLOPT: Nonlinear Optimization [v%s, Built %s, VS%s]\n",PACKAGE_VERSION,__DATE__,vbuf);              
    mexPrintf("  - Released under the GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html\n");
    mexPrintf("  - Source available from: http://ab-initio.mit.edu/wiki/index.php/NLopt\n");

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}
