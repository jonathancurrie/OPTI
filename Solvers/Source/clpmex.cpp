/* CLPMEX - A MATLAB MEX Interface to CLP
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2013
 * www.inverseproblem.co.nz
 */

#include "mex.h"
#include "Coin_C_defines.h"
#include "config_clp_default.h"
#include "config_coinutils_default.h"
#include "ClpSimplex.hpp"
#include "CoinMessageHandler.hpp"
#include "ClpEventHandler.hpp"
#include <exception>
#include <time.h>
#include "opti_util.h"

#ifdef CLP_HAS_ABC
 #include "AbcSimplex.hpp"
 #include <cilk/cilk_api.h>
#endif

using namespace std;

//Argument Enums (in expected order of arguments)
enum {eH, eF, eA, eRL, eRU, eLB, eUB, eOPTS};                   
//PRHS Defines    
#define pH      prhs[eH]
#define pF      prhs[eF]
#define pA      prhs[eA]
#define pRL     prhs[eRL]
#define pRU     prhs[eRU]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pOPTS   prhs[eOPTS]

//Function Prototypes
void printSolverInfo();
void GetIntegerOption(const mxArray *opts, char *name, int *var);
void GetDoubleOption(const mxArray *opts, char *name, double *var);
void checkInputs(const mxArray *prhs[], int nrhs);
void lower(char *str);

//Ctrl-C Detection 
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);

//Message Handler
class DerivedHandler : public CoinMessageHandler {
public:
	virtual int print();  
};
int DerivedHandler::print()
{
	mexPrintf(messageBuffer());
	mexPrintf("\n");
    mexEvalString("drawnow;"); //flush draw buffer
	return 0;
}

//Ctrl-C Event Handler
class DerivedEvent : public ClpEventHandler {
public:
    int event(Event whichEvent)
    {     
        if (utIsInterruptPending()) {
            utSetInterruptPending(false); /* clear Ctrl-C status */
            mexPrintf("\nCtrl-C Detected. Exiting CLP...\n\n");
            return 5; //terminate asap
        }
        else
            return -1; //return ok
    }
    ClpEventHandler * DerivedEvent::clone() const
    {
        return new DerivedEvent(*this);
    }
};

//Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *H, *f, *A, *rl, *ru, *lb = NULL, *ub = NULL;
    
    //Return Args
    double *x, *fval, *exitflag, *iter;
    
    //Options
    double maxtime = 1000, objbias = 0.0;
    double primeTol = 1e-7, dualTol = 1e-7;
    double primeObjLim = COIN_DBL_MAX, dualObjLim = COIN_DBL_MAX;
    int maxiter = 1500, numberRefinements = 0, numThreads = 1, abcState = 0;    
    int printLevel = 0, factorFreq = -999, numPresolvePasses = 5;
    int algorithm = (int)ClpSolve::automatic;
    ClpSolve::PresolveType doPresolve = ClpSolve::presolveOn;
    
    //Internal Vars
    size_t ncon, ndec, i;
    double *sol, *llb, *lub;
    const char *lnames[2] = {"dual_row","dual_col"};
    char errstr[1024], msgstr[128], wrkstr[128];
    mxArray *dRow, *dCol;
        
    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    mwIndex *H_ir, *H_jc;
    mwIndex nnzA = 0;
    int *Hrows = NULL, *Arows = NULL;
    CoinBigIndex *Hcols = NULL, *Acols = NULL;
    
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(CLP_VERSION);
        return;
    }        
    
    //Check Inputs
    checkInputs(prhs,nrhs); 
    
    //Get pointers to Input variables
    H = mxGetPr(pH); H_ir = mxGetIr(pH); H_jc = mxGetJc(pH);
	f = mxGetPr(pF);
	A = mxGetPr(pA); A_ir = mxGetIr(pA); A_jc = mxGetJc(pA);
    rl = mxGetPr(pRL);
    ru = mxGetPr(pRU);
    if(nrhs > eLB && !mxIsEmpty(pLB))
        lb = mxGetPr(pLB); 
    if(nrhs > eUB && !mxIsEmpty(pUB))
        ub = mxGetPr(pUB);
        
    //Get Options if Specified
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        GetIntegerOption(pOPTS, "maxiter", &maxiter);
        GetIntegerOption(pOPTS, "display", &printLevel);               
        GetIntegerOption(pOPTS, "numThreads", &numThreads);
        GetIntegerOption(pOPTS, "abcState", &abcState);
        GetIntegerOption(pOPTS, "numPresolvePasses", &numPresolvePasses);
        GetIntegerOption(pOPTS, "factorFreq", &factorFreq);
        GetIntegerOption(pOPTS, "numberRefinements", &numberRefinements);        
        GetDoubleOption(pOPTS, "maxtime", &maxtime);
        GetDoubleOption(pOPTS, "primalTol", &primeTol);
        GetDoubleOption(pOPTS, "dualTol", &dualTol);
        GetDoubleOption(pOPTS, "primalObjLim", &primeObjLim);
        GetDoubleOption(pOPTS, "dualObjLim", &dualObjLim);
        GetDoubleOption(pOPTS, "objbias", &objbias);
        //Manually process Algorithm
        if(mxGetField(pOPTS,0,"algorithm") && !mxIsEmpty(mxGetField(pOPTS,0,"algorithm"))) {
            if(mxIsChar(mxGetField(pOPTS,0,"algorithm"))) {
                char *str = mxArrayToString(mxGetField(pOPTS,0,"algorithm"));
                lower(str);
                if(!strcmp(str,"dualsimplex"))
                    algorithm = 0;
                else if(!strcmp(str,"primalsimplex"))
                    algorithm = 1;
                else if(!strcmp(str,"primalsimplexorsprint"))
                    algorithm = 2;
                else if(!strcmp(str,"barrier"))
                    algorithm = 3;
                else if(!strcmp(str,"barriernocross"))
                    algorithm = 4;
                else if(!strcmp(str,"automatic"))
                    algorithm = 5;
                else
                    mexErrMsgTxt("Unknown algorithm - options are 'dualsimplex', 'primalsimplex', 'primalsimplexorsprint', 'barrier', 'barriernocross' or 'automatic'");
                mxFree(str);
            }
            else
                algorithm = (int)*mxGetPr(mxGetField(pOPTS,0,"algorithm"));
        }
        //Manually process doPresolve (Due to inverted 0/1 on/off)
        if(mxGetField(pOPTS,0,"doPresolve") && !mxIsEmpty(mxGetField(pOPTS,0,"doPresolve"))) {
            if(*mxGetPr(mxGetField(pOPTS,0,"doPresolve")) == 1)
                doPresolve = ClpSolve::presolveOn;
            else
                doPresolve = ClpSolve::presolveOff;
        }
    }
    
    //Ensure Infinite options are limited to COIN limits
    if(mxIsInf(primeObjLim))
        primeObjLim = COIN_DBL_MAX;
    if(mxIsInf(dualObjLim))
        dualObjLim = COIN_DBL_MAX;
    
    //If numThreads = 0, not specified, set to #cilk workers if Abc available
    if(numThreads == 0) {
        #ifdef CLP_HAS_ABC
         numThreads = __cilkrts_get_nworkers();
        #endif
    }
    //Else if user has specified numThreads and > 1, ensure built with Aboca 
    else if(numThreads > 1) {
        #ifndef CLP_HAS_ABC
         mexWarnMsgTxt("This library has not been compiled with Aboca. Your problem will solve in sequential mode only.");
         numThreads = 1;
        #endif
    }                
         
    //Get sizes
    ndec = mxGetNumberOfElements(pF);
    ncon = mxGetM(pA); 

    //If QP, Aboca not supported
    if(!mxIsEmpty(pH)) {
        mwIndex nzH = H_jc[ndec];
        if(nzH > 0 && numThreads != 1)
        {
            mexWarnMsgTxt("Quadratic Problems are not supported with multiple threads in Clp. Your problem will solve in sequential mode only.");
            numThreads = 1;
        }
    }
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);  
    plhs[4] = mxCreateStructMatrix(1,1,2,lnames);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    iter = mxGetPr(plhs[3]);
    dRow = mxCreateDoubleMatrix(ncon,1, mxREAL);
    dCol = mxCreateDoubleMatrix(ndec,1, mxREAL);
    
    try
    {
        //CLP Objects
        ClpSimplex simplex;
        ClpSolve options;          
        DerivedHandler *mexprinter;   
        DerivedEvent *ctrlCEvent;
        
        //Create bounds
        llb = (double*)mxCalloc(ndec,sizeof(double));
        lub = (double*)mxCalloc(ndec,sizeof(double));
        //If no bounds specified, fill with defaults, otherwise copy
        if(lb == NULL) {            
            for(i=0;i<ndec;i++)
                llb[i] = -COIN_DBL_MAX;
        }
        else {
            for(i=0;i<ndec;i++) {
                if(mxIsInf(lb[i]))
                    llb[i] = -COIN_DBL_MAX;
                else
                    llb[i] = lb[i];
            }
        }
        if(ub == NULL) {
            for(i=0;i<ndec;i++)
                lub[i] = COIN_DBL_MAX;
        }
        else {
            for(i=0;i<ndec;i++) {
                if(mxIsInf(ub[i]))
                    lub[i] = COIN_DBL_MAX;
                else
                    lub[i] = ub[i];
            }
        }
        //Convert Linear Constraints to COIN Format
        if(ncon) {
            //Convert Indicies
            nnzA = A_jc[ndec];
            Arows = (int*)mxCalloc(nnzA,sizeof(int));
            Acols = (CoinBigIndex*)mxCalloc(ndec+1,sizeof(CoinBigIndex));
            //Assign Convert Data Type Vectors
            for(i = 0; i <= ndec; i++)
               Acols[i] = (CoinBigIndex)A_jc[i];
            for(i = 0; i < nnzA; i++)
               Arows[i] = (int)A_ir[i];
        }
        
        //Load Problem Into Selected Solver
        simplex.loadProblem((int)ndec,(int)ncon,Acols,Arows,A,llb,lub,f,rl,ru);
        
        //Add Quadratic Objective if specified
        if(!mxIsEmpty(pH)) {           
            //Convert Indicies
            mwIndex nzH = H_jc[ndec];
            if(nzH) {           
                Hrows = (int*)mxCalloc(nzH,sizeof(int));
                Hcols = (CoinBigIndex*)mxCalloc(ndec+1,sizeof(CoinBigIndex));
                //Assign Convert Data Type Vectors
                for(i = 0; i <= ndec; i++)
                   Hcols[i] = (CoinBigIndex)H_jc[i];
                for(i = 0; i < nzH; i++)
                   Hrows[i] = (int)H_ir[i];
                //Load QuadObj into selected solver
                simplex.loadQuadraticObjective((int)ndec,Hcols,Hrows,H);
            }
        }
        
        if(numThreads > 1) {
            sprintf(msgstr,"Parallel (Aboca) ");
            sprintf(wrkstr," with %d threads",numThreads);
        }
        else {
            sprintf(msgstr,"");
            sprintf(wrkstr,"");
        }
        
        if(printLevel) {
            mexPrintf("\n------------------------------------------------------------------\n");
            mexPrintf(" This is CLP v%s\n Author: John J. Forrest\n\n",CLP_VERSION); 
            switch((ClpSolve::SolveType)algorithm) {
                case ClpSolve::useDual: mexPrintf(" This run uses the %sDual Simplex Solver%s\n\n",msgstr,wrkstr); break;
                case ClpSolve::usePrimal: mexPrintf(" This run uses the %sPrimal Simplex Solver%s\n\n",msgstr,wrkstr); break;
                case ClpSolve::usePrimalorSprint: mexPrintf(" This run uses the %sPrimal or Sprint Simplex Solver%s\n\n",msgstr,wrkstr); break;
                case ClpSolve::useBarrier: mexPrintf(" This run uses the Barrier Interior Point Solver (with Simplex Crossover)%s\n\n",wrkstr); break;
                case ClpSolve::useBarrierNoCross: mexPrintf(" This run uses the Barrier Interior Point Solver (no Crossover)%s\n\n",wrkstr); break;
                case ClpSolve::automatic: mexPrintf(" This run uses an Automatically Chosen %ssolver%s\n\n",msgstr,wrkstr); break;
            }            
            mexPrintf(" Problem Properties:\n # Decision Variables: %6d\n # Linear Constraints: %6d [%d nz]\n",ndec,ncon,nnzA);
            mexPrintf("------------------------------------------------------------------\n");
            mexEvalString("drawnow;");
        }
       
        //Set Options   
        simplex.setPrimalTolerance(primeTol);
        simplex.setDualTolerance(dualTol);
        simplex.setPrimalObjectiveLimit(primeObjLim);
        simplex.setDualObjectiveLimit(dualObjLim);
        simplex.setMaximumIterations(maxiter);
        simplex.setMaximumSeconds(maxtime);
        simplex.setObjectiveOffset(-objbias);
        if(factorFreq = -999) //not set
            factorFreq = 100 + ncon/50; //J.Forrest Heuristic
        simplex.setFactorizationFrequency(factorFreq);
        simplex.setNumberRefinements(numberRefinements);
        if(printLevel) {
            mexprinter = new DerivedHandler();
            mexprinter->setLogLevel(printLevel);
            simplex.passInMessageHandler(mexprinter);
        }
        //Add Event Handler for Ctrl+C
        ctrlCEvent = new DerivedEvent();  
        simplex.passInEventHandler(ctrlCEvent);    

        //Setup ClpSolve Options
        options.setSolveType((ClpSolve::SolveType)algorithm);
        options.setPresolveType((ClpSolve::PresolveType)doPresolve,numPresolvePasses);    

        //If we are solving using Aboca, convert and solve
        #ifdef CLP_HAS_ABC
            if(numThreads > 1) {
                //Set Number of Workers Cilk can Use    
                sprintf(wrkstr,"%d",numThreads);
                __cilkrts_end_cilk(); //ensure cilk rt has stopped 
                int exNW = __cilkrts_get_nworkers();
                __cilkrts_set_param("nworkers", wrkstr); 
                //Convert to AbcSimplex Object
                AbcSimplex parsimplex(simplex);    
                //Set Partition Size
                if(abcState != 0)
                    parsimplex.setAbcState(abcState);
                else
                    parsimplex.setAbcState(min(8,2*numThreads));               
                //Solve Model with Optional Presolve            
                parsimplex.initialSolve(options);        
                //Assign Return Arguments
                sol = parsimplex.primalColumnSolution();
                if(sol != NULL) {
                    memcpy(x,sol,ndec*sizeof(double));
                    *fval = parsimplex.objectiveValue();
                    *exitflag = parsimplex.status();
                    *iter = parsimplex.numberIterations();
                    memcpy(mxGetPr(dRow),parsimplex.dualRowSolution(),ncon*sizeof(double));
                    memcpy(mxGetPr(dCol),parsimplex.dualColumnSolution(),ndec*sizeof(double));
                }
                //Return cilk to default number of workers
                __cilkrts_end_cilk(); //ensure cilk rt has stopped 
                sprintf(wrkstr,"%d",exNW);
                __cilkrts_set_param("nworkers", wrkstr);
            }
            else {
                //Solve Model Sequentially with Optional Presolve            
                simplex.initialSolve(options);        
                //Assign Return Arguments
                sol = simplex.primalColumnSolution();
                if(sol != NULL) {
                    memcpy(x,sol,ndec*sizeof(double));
                    *fval = simplex.objectiveValue();
                    *exitflag = simplex.status();
                    *iter = simplex.numberIterations();
                    memcpy(mxGetPr(dRow),simplex.dualRowSolution(),ncon*sizeof(double));
                    memcpy(mxGetPr(dCol),simplex.dualColumnSolution(),ndec*sizeof(double));
                }
            }
        #else
            //Solve Model Sequentially with Optional Presolve            
            simplex.initialSolve(options);        
            //Assign Return Arguments
            sol = simplex.primalColumnSolution();
            if(sol != NULL) {
                memcpy(x,sol,ndec*sizeof(double));
                *fval = simplex.objectiveValue();
                *exitflag = simplex.status();
                *iter = simplex.numberIterations();
                memcpy(mxGetPr(dRow),simplex.dualRowSolution(),ncon*sizeof(double));
                memcpy(mxGetPr(dCol),simplex.dualColumnSolution(),ndec*sizeof(double));
            }
        #endif
            
        //Assign dual solution
        if(sol != NULL) {
            mxSetField(plhs[4],0,lnames[0],dRow);
            mxSetField(plhs[4],0,lnames[1],dCol);
        }
        
        if(printLevel){            
            switch((int)*exitflag) {
                case -1: mexPrintf("\n *** TERMINATION: PROBABLY INFEASIBLE ***\n"); break;
                case  0: mexPrintf("\n *** SUCCESSFUL TERMINATION: PROVEN OPTIMAL ***\n"); break;
                case  1: mexPrintf("\n *** TERMINATION: PROVEN PRIMAL INFEASIBLE ***\n"); break;
                case  2: mexPrintf("\n *** TERMINATION: PROVEN DUAL INFEASIBLE ***\n"); break;
                case  3: mexPrintf("\n *** MAXIMUM ITERATIONS OR TIME REACHED ***\n"); break;
                case  4: mexPrintf("\n *** TERMINATION: ERRORS PRESENT ***\n"); break;
                case  5: mexPrintf("\n *** USER EXITED ***\n"); break;
                default: mexPrintf("\n *** TERMINATION: UNKNOWN REASON (CODE: %d) ***\n",(int)*exitflag); break;
            }            
            if(sol != NULL)
                mexPrintf(" Final Objective Value: %9.6g\n In %3.0f iterations\n",*fval,*iter);
            
            mexPrintf("------------------------------------------------------------------\n\n");
        }
       
        //Clean up memory  
        if(printLevel)
           delete mexprinter;
        mxFree(llb); mxFree(lub);        
        if(Hrows) {mxFree(Hrows); Hrows = NULL;}
        if(Hcols) {mxFree(Hcols); Hcols = NULL;}
        if(Arows) {mxFree(Arows); Arows = NULL;}
        if(Acols) {mxFree(Acols); Acols = NULL;}
        
    }
    //Error Handling
    catch(CoinError e) {
        sprintf(errstr,"Caught Coin Error: %s",e.message());
        mexErrMsgTxt(errstr);
    }
    catch(exception& e) {
        sprintf(errstr,"Caught CLP Error: %s",e.what());           
        mexErrMsgTxt(errstr);
    }  
    catch(...) {
        mexErrMsgTxt("Unknown Error Running CLP");
    }
}               


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, ncon;
    
    //Correct number of inputs
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 5 arguments to clp (H, f, A, rl, ru)"); 
    
    //Check we have an objective
    if(mxIsEmpty(pH) && mxIsEmpty(pF))
        mexErrMsgTxt("You must supply an objective function!");
    if(mxIsEmpty(pF))
        mexErrMsgTxt("You must supply f (linear objective vector)!");
    
    //Check we have some constraints
    if(nrhs <= eLB) {
        if(mxIsEmpty(pA))
            mexErrMsgTxt("You have not supplied any constraints!");
    }
    else {
        if(nrhs > eUB) {
            if(mxIsEmpty(pA) && mxIsEmpty(pLB) && mxIsEmpty(pUB))
                mexErrMsgTxt("You have not supplied any constraints!");
        }
        else {
            if(mxIsEmpty(pA) && mxIsEmpty(pLB))
                mexErrMsgTxt("You have not supplied any constraints!");
        }
    }
   
    //Check options is a structure
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && !mxIsStruct(pOPTS))
        mexErrMsgTxt("The options argument must be a structure!");
    
    //Get Sizes
    ndec = mxGetNumberOfElements(pF);
    ncon = mxGetM(pA);
    
    //Check Constraint Pairs
    if(ncon && mxIsEmpty(pRL))
        mexErrMsgTxt("When A is specified, rl must not be empty!");
    if(ncon && mxIsEmpty(pRU))
        mexErrMsgTxt("When A is specified, ru must not be empty!");
    
    //Check Data Types and Sparsity
    if(!mxIsEmpty(pH) && (!mxIsSparse(pH) || !mxIsDouble(pH) || mxIsComplex(pH)))
        mexErrMsgTxt("H must be a real, sparse, double matrix");
    if(mxIsSparse(pF) || !mxIsDouble(pF) || mxIsComplex(pF))
        mexErrMsgTxt("f must be a real, dense, double vector");
    if(!mxIsEmpty(pA) && (!mxIsSparse(pA) || !mxIsDouble(pA) || mxIsComplex(pA)))
        mexErrMsgTxt("A must be a real, sparse, double matrix");
    if(!mxIsEmpty(pRL) && (mxIsSparse(pRL) || !mxIsDouble(pRL) || mxIsComplex(pRL)))
        mexErrMsgTxt("rl must be a real, dense, double vector");
    if(!mxIsEmpty(pRU) && (mxIsSparse(pRU) || !mxIsDouble(pRU) || mxIsComplex(pRU)))
        mexErrMsgTxt("ru must be a real, dense, double vector");
    if(nrhs > eLB && !mxIsEmpty(pLB) && (mxIsSparse(pLB) || !mxIsDouble(pLB) || mxIsComplex(pLB)))
        mexErrMsgTxt("lb must be a real, dense, double vector");
    if(nrhs > eUB && !mxIsEmpty(pUB) && (mxIsSparse(pUB) || !mxIsDouble(pUB) || mxIsComplex(pUB)))
        mexErrMsgTxt("ub must be a real, dense, double vector");
    
    //Check Sizes
    if(!mxIsEmpty(pH)) {
        if(mxGetN(pH) != ndec || mxGetM(pH) != ndec)
            mexErrMsgTxt("H has incompatible dimensions");
    }
    if(ncon) {
        if(mxGetN(pA) != ndec)
            mexErrMsgTxt("A has incompatible dimensions");
        if(mxGetNumberOfElements(pRL) != ncon)
            mexErrMsgTxt("rl has incompatible dimensions");
        if(mxGetNumberOfElements(pRU) != ncon)
            mexErrMsgTxt("ru has incompatible dimensions");
    }
    if(nrhs > eLB && !mxIsEmpty(pLB) && (mxGetNumberOfElements(pLB) != ndec))
        mexErrMsgTxt("lb has incompatible dimensions");
    if(nrhs > eUB && !mxIsEmpty(pUB) && (mxGetNumberOfElements(pUB) != ndec))
        mexErrMsgTxt("ub has incompatible dimensions");   
    
    //Ensure H is lower triangular
    if(!mxIsEmpty(pH)) {
        mwIndex *jc = mxGetJc(pH);
        mwIndex *ir = mxGetIr(pH);
        mwIndex k = 0;
        mwSize n = mxGetN(pH);
        for(mwIndex i = 0; i < n; i++) {
            mwIndex start = jc[i];
            mwIndex stop = jc[i+1];
            for(mwIndex j = start; j < stop; j++) {
                if(i > ir[k++])
                    mexErrMsgTxt("H is not symmetric lower triangular");
            }
        }        
    } 
}

//Print Solver Information
void printSolverInfo()
{     
    char vbuf[6];   
    mexPrintf("\n-----------------------------------------------------------\n");
    #ifdef CLP_HAS_ABC
        mexPrintf(" CLP: COIN-OR Linear Programming [v%s, Built %s] with ABOCA (A Bit of Coin-Or Accelerated)\n",CLP_VERSION,__DATE__);
    #else
        getVSVer(vbuf);  
        mexPrintf(" CLP: COIN-OR Linear Programming [v%s, Built %s, VS%s]\n",CLP_VERSION,__DATE__,vbuf);
    #endif
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/Clp\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - CoinUtils [v%s] (Eclipse Public License)\n",COINUTILS_VERSION);

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
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

//Convert input string to lowercase
void lower(char *str)
{
    int i = 0;
    while(str[i]) {
        str[i] = tolower(str[i]);  
        i++;
    }
}