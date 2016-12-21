/* CBCMEX - A MATLAB MEX Interface to CBC
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2013
 * www.inverseproblem.co.nz
 */

#include "mex.h"
#include "Coin_C_defines.h"
#include "config_clp_default.h"
#include "config_cbc_default.h"
#include "config_cgl_default.h"
#include "config_osi_default.h"
#include "config_coinutils_default.h"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcSOS.hpp"
#include "CoinMessageHandler.hpp"
#include "CbcEventHandler.hpp"
#include "CbcOrClpParam.hpp"
#include "opti_util.h"

#include <exception>
#include <string>
#include <list>

using namespace std;

//Argument Enums (in expected order of arguments)
enum {eH, eF, eA, eRL, eRU, eLB, eUB, eXTYPE, eSOS, eX0, eOPTS};                   
//PRHS Defines    
#define pH      prhs[eH]
#define pF      prhs[eF]
#define pA      prhs[eA]
#define pRL     prhs[eRL]
#define pRU     prhs[eRU]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pXTYPE  prhs[eXTYPE]
#define pSOS    prhs[eSOS]
#define pX0     prhs[eX0]
#define pOPTS   prhs[eOPTS]

//Function Prototypes
void printSolverInfo();
bool isCharOption(const mxArray *opts, char *name);
bool isRealOption(const mxArray *opts, char *name);
void GetIntegerOption(const mxArray *opts, char *name, int *var);
void GetDoubleOption(const mxArray *opts, char *name, double *var);
int getIntegerOption(const mxArray *opts, char *name);
double getDoubleOption(const mxArray *opts, char *name);
void setupParameterList(const mxArray *opts, std::list<std::string>& par_list);
void checkInputs(const mxArray *prhs[], int nrhs);

//Ctrl-C Detection
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);

//Message Handler
class DerivedHandler : public CoinMessageHandler {
public:
    DerivedHandler() : CoinMessageHandler() {inSmall_ = false;}
	virtual int print() ;
protected:
    bool inSmall_; 
};
int DerivedHandler::print()
{
    if (!inSmall_) {
        mexPrintf(messageBuffer());
        mexPrintf("\n");
        mexEvalString("drawnow;"); //flush draw buffer

        if (currentMessage().externalNumber()==38)
            inSmall_=true;
    }

    if (inSmall_) {
        if (currentMessage().externalNumber()==35 || currentMessage().externalNumber()==41)
            inSmall_=false;
    }
    return 0;
}


//Ctrl-C Event Handler
class DerivedEvent : public CbcEventHandler {
public:
     virtual CbcAction event(CbcEvent whichEvent);
     virtual CbcEventHandler * clone() const ;
};
CbcEventHandler::CbcAction DerivedEvent::event(CbcEvent whichEvent)
{
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting CBC...\n\n");
        return stop; //terminate asap
    }
    else
        return noAction; //return ok
}
CbcEventHandler * DerivedEvent::clone() const
{
     return new DerivedEvent(*this);
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *H = NULL, *f, *A = NULL, *rl = NULL, *ru = NULL, *lb = NULL, *ub = NULL, *x0 = NULL;
    char *xtype = NULL;
    
    //Return Args
    double *x, *fval, *exitflag, *nodes, *contobj;
    
    //Options
    int maxnodes = 100000, printLevel = 0;  
    double maxtime = 1000, intTol = 1e-5, primalTol = 1e-7, dualTol = 1e-7;
    double gap = 0, fracgap = 0;
    double objbias = 0.0;
    
    //Internal Vars
    size_t ncon = 0, ndec = 0, nint = 0, nbin = 0, i, j;
    double *sol, *llb, *lub, *lrl = NULL, *lru = NULL;
    int iters = 0;
    char errstr[1024];
    
    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    mwIndex *H_ir, *H_jc;
    mwIndex nnzA = 0;
    int *Hrows = NULL, *Arows = NULL;
    CoinBigIndex *Hcols = NULL, *Acols = NULL;
    
    //SOS
    CbcObject ** objects;
    int no_sets = 0;
    char *sostype;
    double *sosind, *soswt;
    int no_entries, type, *isosind; 
    
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(CBC_VERSION);
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
    if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE))
        xtype = mxArrayToString(pXTYPE);
    if(nrhs > eX0 && !mxIsEmpty(pX0))
        x0 = mxGetPr(pX0);

    //Get sizes
    ndec = mxGetNumberOfElements(pF);
    ncon = mxGetM(pA);   
    
    //Get Easy Options if Specified
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        GetIntegerOption(pOPTS, "display", &printLevel);        
        GetIntegerOption(pOPTS, "maxnodes", &maxnodes);
        GetDoubleOption(pOPTS, "maxtime", &maxtime);
        GetDoubleOption(pOPTS, "primalTol", &primalTol);
        GetDoubleOption(pOPTS, "dualTol", &dualTol);
        GetDoubleOption(pOPTS, "intTol", &intTol);
        GetDoubleOption(pOPTS, "allowableGap", &gap);
        GetDoubleOption(pOPTS, "allowableFracGap", &fracgap);      
        GetDoubleOption(pOPTS, "objbias", &objbias);  
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
    nodes = mxGetPr(plhs[3]);
    contobj = mxGetPr(plhs[4]);
    
    try
    {
        //Objects
        OsiClpSolverInterface OSImodel;        
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
        //Modify bounds based on binary constraints
        if(xtype) {
            for(i = 0; i < ndec; i++) {
                switch(tolower(xtype[i]))
                {
                    case 'b':
                        //Enforce binary bounds if not specified
                        if(llb[i] == -COIN_DBL_MAX)
                            llb[i] = 0;
                        if(lub[i] == COIN_DBL_MAX)
                            lub[i] = 1;
                        break;
                }
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
            
            //Copy and process row bounds
            lrl = (double*)mxCalloc(ncon,sizeof(double));
            lru = (double*)mxCalloc(ncon,sizeof(double)); 
            for(i = 0; i < ncon; i++) {
                if(mxIsInf(rl[i]))
                    lrl[i] = -COIN_DBL_MAX;
                else
                    lrl[i] = rl[i];
                
                if(mxIsInf(ru[i]))
                    lru[i] = COIN_DBL_MAX;
                else
                    lru[i] = ru[i];
            }
        }
        
        //Load Problem
        OSImodel.setDblParam(OsiObjOffset, objbias); 
        OSImodel.loadProblem((int)ndec,(int)ncon,Acols,Arows,A,llb,lub,f,lrl,lru);

        //Add Integer and Binary Constraints
        if(xtype) {
            for(i = 0; i < ndec; i++) {
                switch(tolower(xtype[i]))
                {
                    case 'c': break;
                    case 'i': 
                        nint++;
                        OSImodel.setInteger((int)i); 
                        break;
                    case 'b':
                        nbin++;
                        OSImodel.setInteger((int)i); 
                        break;
                    default:
                        throw exception("Unknown xtype, only 'C', 'I' and 'B' are accepted");
                }
            } 
        }  
        
        //Add Quadratic Objective if specified  NOT WORKING
        if(!mxIsEmpty(pH)) {    
            mexWarnMsgTxt("Quadratic Objectives Are Not Currently Supported in this Interface for Cbc! Problem will be treated as a (MI)LP only.");
            /*//Convert Indicies
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
                OSImodel.getModelPtr()->loadQuadraticObjective((int)ndec,Hcols,Hrows,H);
            }*/
        }
        
        //Load Problem into CBC
        CbcModel model(OSImodel);
        //Setup SOS
        if(nrhs > eSOS && !mxIsEmpty(pSOS)) {
            no_sets = (int)mxGetNumberOfElements(mxGetField(pSOS,0,"type"));
            if(no_sets > 0) {
                //Collect Types
                sostype = mxArrayToString(mxGetField(pSOS,0,"type"));
                //Allocate Set Memory
                objects = new CbcObject * [no_sets];
                //Copy in SOS data, creating CbcSOS objects as we go
                for(i=0;i<no_sets;i++) {
                    type = (int)sostype[i] - '0'; //convert to numerical representation
                    if(type < 1 || type > 2)
                        throw exception("Only SOS of type '1' and '2' are supported");
                    no_entries = (int)mxGetNumberOfElements(mxGetCell(mxGetField(pSOS,0,"index"),i));
                    //Create sosind memory and copy in
                    isosind = new int[no_entries];
                    sosind = mxGetPr(mxGetCell(mxGetField(pSOS,0,"index"),i));
                    for(j=0;j<no_entries;j++)
                        isosind[j] = (int)sosind[j]-1;                
                    //Get soswt
                    soswt = mxGetPr(mxGetCell(mxGetField(pSOS,0,"weight"),i));
                    //Create Set object
                    objects[i] = new CbcSOS(&model,no_entries,isosind,soswt,(int)i,type);
                    delete []isosind; //free memory as we go round, CoinSet copies internally
                }

                //Add objects to model
                model.addObjects(no_sets,objects); 
                //Priorities??

                //Delete objects
                for(i=0;i<no_sets;i++)
                    delete objects[i];
                delete [] objects;
                //Delete type string
                mxFree(sostype);
            }
        }
        //Initial Guess If Specified
        if(x0)
            model.setHotstartSolution(x0);
        
        if(printLevel) {
            mexPrintf("\n------------------------------------------------------------------\n");
            mexPrintf(" This is CBC v%s\n Author: John J. Forrest\n\n",CBC_VERSION);            
            mexPrintf(" Problem Properties:\n # Decision Variables: %6d [%d Integer, %d Binary]\n # Linear Constraints: %6d [%d nz]\n",ndec,nint,nbin,ncon,nnzA);
            mexPrintf("------------------------------------------------------------------\n");
            mexEvalString("drawnow;");
        }
        
        //Initialize Options
        CbcMain0(model);
        
        //Set CbcModel Options        
        model.setMaximumNodes(maxnodes);
		model.setMaximumSeconds(maxtime);
        model.solver()->setDblParam(OsiPrimalTolerance, primalTol);
        model.solver()->setDblParam(OsiPrimalTolerance, dualTol);
        model.setIntegerTolerance(intTol);        
        model.setAllowableGap(gap);
        model.setAllowableFractionGap(fracgap);
        //Message Handler
        if(printLevel) {
            mexprinter = new DerivedHandler();            
            model.passInMessageHandler(mexprinter);
            model.messageHandler()->setLogLevel(printLevel);
            model.messageHandler()->setLogLevel(0,printLevel);
        }
        else
            model.messageHandler()->setLogLevel(0);

        //Add Event Handler for Ctrl+C
        ctrlCEvent = new DerivedEvent();      
        model.passInEventHandler(ctrlCEvent); 
        
        //Set Remaining Command Line Options [based on GamsCoinCbc.cpp - GAMSLINKS]
        if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
            std::list<std::string> par_list; 
            setupParameterList(pOPTS, par_list); 
            size_t par_list_length=par_list.size(); 
            const char** cbc_args=new const char*[par_list_length+2]; 
            cbc_args[0]="OPTI"; i=1; 
            for (std::list<std::string>::iterator it(par_list.begin()); it!=par_list.end(); ++it, ++i) 
                cbc_args[i]=it->c_str(); 
            cbc_args[i++]="-quit"; 

            //Solve using CBC + CLP
            CbcMain1((int)(par_list_length+2),cbc_args,model);
            delete[] cbc_args;
        }
        else {
            const char *cbc_args[] = {"OPTI","-solve","-quit"};
            CbcMain1(3,cbc_args,model);
        }
        
        //Assign Return Arguments
        sol = (double*)model.getColSolution();
        
        if(sol != NULL) {
            memcpy(x,sol,ndec*sizeof(double));
            *fval = model.getObjValue();
            //Check for proven infeasible, otherwise use secondary status
            if(model.status() == 0 && model.isProvenInfeasible())
                *exitflag = 8.;
            else
                *exitflag = (double)model.secondaryStatus();            
            *nodes = model.getNodeCount();
            *contobj = model.getContinuousObjective();
            iters = model.getIterationCount();
        }
        
        if(printLevel){            
            switch((int)*exitflag) {
                case  0: mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n"); break;
                case  1: mexPrintf("\n *** TERMINATION: LINEAR RELAXATION INFEASIBLE ***\n"); break;
                case  2: mexPrintf("\n *** TERMINATION: GAP REACHED ***\n"); break;
                case  3: mexPrintf("\n *** MAXIMUM NODES REACHED ***\n"); break;
                case  4: mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;                
                case  5: mexPrintf("\n *** USER EXITED ***\n"); break;
                case  6: mexPrintf("\n *** TERMINATION: #SOLUTIONS REACHED ***\n"); break;
                case  7: mexPrintf("\n *** TERMINATION: LINEAR RELAXATION UNBOUNDED ***\n"); break;
                case  8: mexPrintf("\n *** TERMINATION: PROVEN INFEASIBLE ***\n"); break;
                default: mexPrintf("\n *** TERMINATION: UNKNOWN REASON (CODE: %d) ***\n",(int)*exitflag); break;
            }     
            double calcrgap = abs(*contobj-*fval)/(1e-1 + abs(*fval));
            if(sol != NULL) {
                mexPrintf(" Objective Value: %1.6g\n",*fval);
                if(calcrgap < 1e300)
                    mexPrintf(" Gap:             %1.6g (%1.2f%%)\n",abs(*contobj-*fval),calcrgap*100);
                mexPrintf(" Searched:        %1.0f nodes [%1d LP iterations]\n",*nodes,iters);
            }            
            mexPrintf("------------------------------------------------------------------\n\n");
        }

        //Clean up memory
        mxFree(llb); mxFree(lub);
        if(lrl)   mxFree(lrl);   lrl = NULL;
        if(lru)   mxFree(lru);   lru = NULL;
        if(xtype) mxFree(xtype); xtype = NULL;
        if(printLevel)
            delete mexprinter;
    }
    //Error Handling
    catch(CoinError e) {
        sprintf(errstr,"Caught Coin Error: %s",e.message());
        mexErrMsgTxt(errstr);
    }
    catch(exception& e) {
        sprintf(errstr,"Caught CBC Error: %s",e.what());           
        mexErrMsgTxt(errstr);
    }  
    catch(...) {
        mexErrMsgTxt("Unknown Error Running CBC");
    } 
}               


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, ncon;
    
    //Correct number of inputs
    if(nrhs < 4)
        mexErrMsgTxt("You must supply at least 5 arguments to cbc (H, f, A, rl, ru)"); 
    
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
    if(nrhs > eX0 && !mxIsEmpty(pX0) && (mxIsSparse(pX0) || !mxIsDouble(pX0) || mxIsComplex(pX0)))
        mexErrMsgTxt("x0 must be a real, dense, double vector");
    
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
    if(nrhs > eX0 && !mxIsEmpty(pX0) && (mxGetNumberOfElements(pX0) != ndec))
        mexErrMsgTxt("x0 has incompatible dimensions");
    
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
    
     //Check SOS structure
    if(nrhs > eSOS && !mxIsEmpty(pSOS)) {
        if(!mxIsStruct(pSOS))
            mexErrMsgTxt("The SOS argument must be a structure!");       
        if(mxGetFieldNumber(pSOS,"type") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'type'");
        if(mxGetFieldNumber(pSOS,"index") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'index'");
        if(mxGetFieldNumber(pSOS,"weight") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'weight'");
        //Ensure type is char array
        if(!mxIsChar(mxGetField(pSOS,0,"type")))
            mexErrMsgTxt("sos.type should be a char array");
        //Check multiple sets length
        int no_sets = (int)mxGetNumberOfElements(mxGetField(pSOS,0,"type")); 
        if(no_sets > 1) {
            if(!mxIsCell(mxGetField(pSOS,0,"index")) || mxIsEmpty(mxGetField(pSOS,0,"index")))
                mexErrMsgTxt("sos.index must be a cell array, and not empty!");
            if(!mxIsCell(mxGetField(pSOS,0,"weight")) || mxIsEmpty(mxGetField(pSOS,0,"weight")))
                mexErrMsgTxt("sos.weight must be a cell array, and not empty!");
            if(mxGetNumberOfElements(mxGetField(pSOS,0,"index")) != no_sets)
                mexErrMsgTxt("sos.index cell array is not the same length as sos.type!");
            if(mxGetNumberOfElements(mxGetField(pSOS,0,"weight")) != no_sets)
                mexErrMsgTxt("sos.weight cell array is not the same length as sos.type!");        
        }
    }
    
    //Check xtype
    if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE) && (!mxIsChar(pXTYPE) || mxIsEmpty(pXTYPE) || mxGetNumberOfElements(pXTYPE) < ndec))
        mexErrMsgTxt("xtype should be a char array with ndec elements");   
}

//Print Solver Information
void printSolverInfo()
{    
    char vbuf[6]; getVSVer(vbuf);  
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" CBC: COIN-OR Branch and Cut [v%s, Built %s, VS%s]\n",CBC_VERSION,__DATE__,vbuf);
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/Cbc\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - CGL    [v%s] (Eclipse Public License)\n",CGL_VERSION);
    mexPrintf("  - CLP    [v%s] (Eclipse Public License)\n",CLP_VERSION);
    mexPrintf("  - CoinUtils [v%s] (Eclipse Public License)\n",COINUTILS_VERSION);
    mexPrintf("  - OSI    [v%s] (Eclipse Public License)\n",OSI_VERSION);
    
    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}

//Option Getting Methods
bool isCharOption(const mxArray *opts, char *name)
{
    return mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)) && mxIsChar(mxGetField(opts,0,name));
}
bool isRealOption(const mxArray *opts, char *name)
{
    return mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)) && mxIsDouble(mxGetField(opts,0,name));
}
void GetIntegerOption(const mxArray *opts, char *name, int *var)
{
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
        *var = (int)*mxGetPr(mxGetField(opts,0,name));
}
int getIntegerOption(const mxArray *opts, char *name)
{
    return (int)*mxGetPr(mxGetField(opts,0,name));
}
void GetDoubleOption(const mxArray *opts, char *name, double *var)
{
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
        *var = *mxGetPr(mxGetField(opts,0,name));
}
double getDoubleOption(const mxArray *opts, char *name)
{
    return *mxGetPr(mxGetField(opts,0,name));
}

//Parameter Command Line Arguments [From GamsCoinCbc.cpp]
void setupParameterList(const mxArray *opts, std::list<std::string>& par_list) 
{ 
    char buffer[255];
        
	// LP parameters	
	if (isRealOption(opts,"idiotCrash")) {
		par_list.push_back("-idiotCrash");
		sprintf(buffer, "%d", getIntegerOption(opts,"idiotCrash"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"sprintCrash")) {
		par_list.push_back("-sprintCrash");
		sprintf(buffer, "%d", getIntegerOption(opts,"sprintCrash"));
		par_list.push_back(buffer);
	} 
	if (isCharOption(opts,"crash")) {
        char* value = mxArrayToString(mxGetField(opts,0,"crash"));
		par_list.push_back("-crash");
		par_list.push_back(value);
        mxFree(value);				
	}
    if (isCharOption(opts,"factorization")) {
        char* value = mxArrayToString(mxGetField(opts,0,"factorization"));
		par_list.push_back("-factorization");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isRealOption(opts,"maxFactor")) {
		par_list.push_back("-maxFactor");
		sprintf(buffer, "%d", getIntegerOption(opts,"maxFactor"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"crossover")) { // should be revised if we can do quadratic
        char* value = mxArrayToString(mxGetField(opts,0,"crossover"));
		par_list.push_back("-crossover");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"dualPivot")) {
        char* value = mxArrayToString(mxGetField(opts,0,"dualPivot"));
		par_list.push_back("-dualPivot");
		par_list.push_back(value);
        mxFree(value);       			
	}
	if (isCharOption(opts,"primalPivot")) {
        char* value = mxArrayToString(mxGetField(opts,0,"primalPivot"));
		par_list.push_back("-primalPivot");
		par_list.push_back(value);
        mxFree(value);        			
	}	
	if (isCharOption(opts,"perturbation")) {
        char* value = mxArrayToString(mxGetField(opts,0,"perturbation"));
		par_list.push_back("-perturbation");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"scaling")) {
        char* value = mxArrayToString(mxGetField(opts,0,"scaling"));
		par_list.push_back("-scaling");
		par_list.push_back(value);
        mxFree(value);
	}				
	if (isCharOption(opts,"presolve")) {
        char* value = mxArrayToString(mxGetField(opts,0,"presolve"));
		par_list.push_back("-presolve");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isRealOption(opts,"preTolerance")) {
		par_list.push_back("-preTolerance");
		sprintf(buffer, "%g", getDoubleOption(opts,"preTolerance"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"passPresolve")) {
		par_list.push_back("-passPresolve");
		sprintf(buffer, "%d", getIntegerOption(opts,"passPresolve"));
		par_list.push_back(buffer);
	}
    if (isRealOption(opts,"primalWeight")) {
		par_list.push_back("-primalWeight");
		sprintf(buffer, "%g", getDoubleOption(opts,"primalWeight"));
		par_list.push_back(buffer);
	}
    if (isRealOption(opts,"dualBound")) {
		par_list.push_back("-dualBound");
		sprintf(buffer, "%g", getDoubleOption(opts,"dualBound"));
		par_list.push_back(buffer);
	}

	// MIP parameters	
	if (isRealOption(opts,"strategy")) {
		par_list.push_back("-strategy");
		sprintf(buffer, "%d", getIntegerOption(opts,"strategy"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"maxSolutions")) {
		par_list.push_back("-maxSolutions");
		sprintf(buffer, "%d", getIntegerOption(opts,"maxSolutions"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"strongBranching")) {
		par_list.push_back("-strongBranching");
		sprintf(buffer, "%d", getIntegerOption(opts,"strongBranching"));
		par_list.push_back(buffer);
	}			
	if (isRealOption(opts,"trustPseudoCosts")) {
		par_list.push_back("-trustPseudoCosts");
		sprintf(buffer, "%d", getIntegerOption(opts,"trustPseudoCosts"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"cutDepth")) {
		par_list.push_back("-cutDepth");
		sprintf(buffer, "%d", getIntegerOption(opts,"cutDepth"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"maxCutPassesRoot")) {
		par_list.push_back("-passCuts");
		sprintf(buffer, "%d", getIntegerOption(opts,"maxCutPassesRoot"));
		par_list.push_back(buffer);
	}
	if (isRealOption(opts,"maxCutPasses")) {
		par_list.push_back("-passTree");
		sprintf(buffer, "%d", getIntegerOption(opts,"maxCutPasses"));
		par_list.push_back(buffer);
	}
	if (isCharOption(opts,"slowCutPasses")) {
		par_list.push_back("-slowCutPasses");
		sprintf(buffer, "%d", getIntegerOption(opts,"slowCutPasses"));
		par_list.push_back(buffer);
	}
    
    //Cut Settings         
    if (isCharOption(opts,"cliqueCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"cliqueCuts"));
		par_list.push_back("-cliqueCuts");
		par_list.push_back(value);
        mxFree(value);
	}  
    if (isCharOption(opts,"flowCoverCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"flowCoverCuts"));
		par_list.push_back("-flowCoverCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"GMICuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"GMICuts"));
		par_list.push_back("-GMICuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"lagomoryCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"lagomoryCuts"));
		par_list.push_back("-lagomoryCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"gomoryCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"gomoryCuts"));
		par_list.push_back("-gomoryCuts");
		par_list.push_back(value);
        mxFree(value);
	}    
    if (isCharOption(opts,"knapsackCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"knapsackCuts"));
		par_list.push_back("-knapsackCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"liftAndProjectCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"liftAndProjectCuts"));
		par_list.push_back("-liftAndProjectCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"mixedIntegerRoundingCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"mixedIntegerRoundingCuts"));
		par_list.push_back("-mixedIntegerRoundingCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"probingCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"probingCuts"));
		par_list.push_back("-probingCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"reduceAndSplitCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"reduceAndSplitCuts"));
		par_list.push_back("-reduceAndSplitCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"residualCapacityCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"residualCapacityCuts"));
		par_list.push_back("-residualCapacityCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"twoMirCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"twoMirCuts"));
		par_list.push_back("-twoMirCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"latwoMirCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"latwoMirCuts"));
		par_list.push_back("-latwoMirCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"zeroHalfCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"zeroHalfCuts"));
		par_list.push_back("-zeroHalfCuts");
		par_list.push_back(value);
        mxFree(value);
	}
    if (isCharOption(opts,"useCuts")) {
		char* value = mxArrayToString(mxGetField(opts,0,"useCuts"));
		par_list.push_back("-cutsOnOff");
		par_list.push_back(value);
        mxFree(value);
	}
    
    //MIP/Cut Settings On/Off
    if (isCharOption(opts,"heuristics")) {
        char* value = mxArrayToString(mxGetField(opts,0,"heuristics"));
		par_list.push_back("-heuristicsOnOff");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"combineSolution")) {
        char* value = mxArrayToString(mxGetField(opts,0,"combineSolution"));
		par_list.push_back("-combineSolution");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"dins")) {
        char* value = mxArrayToString(mxGetField(opts,0,"dins"));
		par_list.push_back("-Dins");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"divingSome")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingSome"));
		par_list.push_back("-DivingSome");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingCoefficient")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingCoefficient"));
		par_list.push_back("-DivingCoefficient");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingFractional")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingFractional"));
		par_list.push_back("-DivingFractional");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingGuided")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingGuided"));
		par_list.push_back("-DivingGuided");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingLineSearch")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingLineSearch"));
		par_list.push_back("-DivingLineSearch");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingPseudoCost")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingPseudoCost"));
		par_list.push_back("-DivingPseudoCost");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"divingVectorLength")) {
        char* value = mxArrayToString(mxGetField(opts,0,"divingVectorLength"));
		par_list.push_back("-DivingVectorLength");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"feasibilityPump")) {
        char* value = mxArrayToString(mxGetField(opts,0,"feasibilityPump"));
		par_list.push_back("-feasibilityPump");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"localTreeSearch")) {
        char* value = mxArrayToString(mxGetField(opts,0,"localTreeSearch"));
		par_list.push_back("-localTreeSearch");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"naiveHeuristics")) {
        char* value = mxArrayToString(mxGetField(opts,0,"naiveHeuristics"));
		par_list.push_back("-naiveHeuristics");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"pivotAndFix")) {
        char* value = mxArrayToString(mxGetField(opts,0,"pivotAndFix"));
		par_list.push_back("-pivotAndFix");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"randomizedRounding")) {
        char* value = mxArrayToString(mxGetField(opts,0,"randomizedRounding"));
		par_list.push_back("-randomizedRounding");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"rens")) {
        char* value = mxArrayToString(mxGetField(opts,0,"rens"));
		par_list.push_back("-Rens");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"rins")) {
        char* value = mxArrayToString(mxGetField(opts,0,"rins"));
		par_list.push_back("-Rins");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"roundingHeuristic")) {
        char* value = mxArrayToString(mxGetField(opts,0,"roundingHeuristic"));
		par_list.push_back("-roundingHeuristic");
		par_list.push_back(value);
        mxFree(value);
	}
    
    //MIP/Cuts Other Options
	if (isRealOption(opts,"vubHeuristic")) {
		par_list.push_back("-vubheuristic");
		sprintf(buffer, "%d", getIntegerOption(opts,"vubHeuristic"));
		par_list.push_back(buffer);
	}
    if (isRealOption(opts,"feasibilityPumpPasses")) {        
		par_list.push_back("-passFeasibilityPump");
		sprintf(buffer, "%d", getIntegerOption(opts,"feasibilityPumpPasses"));
		par_list.push_back(buffer);
	}	    
    if (isRealOption(opts,"increment")) {
		par_list.push_back("-increment");
		sprintf(buffer, "%g", getDoubleOption(opts,"increment"));
		par_list.push_back(buffer);
	}
    if (isCharOption(opts,"greedyHeuristic")) {        
        char* value = mxArrayToString(mxGetField(opts,0,"greedyHeuristic"));
		par_list.push_back("-greedyHeuristic");
		par_list.push_back(value);
        mxFree(value);
	}    
	if (isCharOption(opts,"costStrategy")) {
        char* value = mxArrayToString(mxGetField(opts,0,"costStrategy"));
		par_list.push_back("-costStrategy");
		par_list.push_back(value);
        mxFree(value);
	}	
	if (isCharOption(opts,"nodeStrategy")) {
        char* value = mxArrayToString(mxGetField(opts,0,"nodeStrategy"));
		par_list.push_back("-nodeStrategy");
		par_list.push_back(value);
        mxFree(value);
	}
	if (isCharOption(opts,"preprocess")) {
        char* value = mxArrayToString(mxGetField(opts,0,"preprocess"));
		par_list.push_back("-preprocess");
		par_list.push_back(value);
        mxFree(value);
	} 
	//Algorithm for root node and solve command 
	if (isCharOption(opts,"startalg")) {
        char* value = mxArrayToString(mxGetField(opts,0,"startalg"));
		if (strcmp(value, "primal")==0) {
			par_list.push_back("-primalSimplex");
		} else if (strcmp(value, "dual")==0) {
			par_list.push_back("-dualSimplex");
		} else if (strcmp(value, "barrier")==0) {
			par_list.push_back("-barrier");
		} 
		par_list.push_back("-solve");
	} else
		par_list.push_back("-solve"); 
}
