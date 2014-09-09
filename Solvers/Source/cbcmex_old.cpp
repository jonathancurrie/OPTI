/* CBCMEX - A MATLAB MEX Interface to CBC
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2013
 * www.i2c2.aut.ac.nz
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
#include "CoinModel.hpp"
#include "CoinMessageHandler.hpp"
#include "CbcEventHandler.hpp"
#include "CglProbing.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcSOS.hpp"

#include <exception>


using namespace std;

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
double getStatus(int stat);

//Ctrl-C Detection (Undocumented - Found in gurobi_mex.c!)
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif

//Message Handler
class DerivedHandler : public CoinMessageHandler {
public:
	virtual int print() ;
};
int DerivedHandler::print()
{
	mexPrintf(messageBuffer());
	mexPrintf("\n");
    mexEvalString("drawnow;"); //flush draw buffer
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
    double *f, *A, *rl, *ru, *lb, *ub, *H = NULL;
    int *ixtype;
    char *xtype = NULL;
    
    //Return Args
    double *x, *fval, *exitflag, *iter, *contobj;
    
    //Options
    int maxnodes = 10000, printLevel = 0, debug = 0;  
    double maxtime = 1000, intTol = 1e-5;
    double objc = 0.0;
    
    //Internal Vars
    size_t ncon, ndec;
    size_t i, j, k;
    const double *sol;
    double *llb, *lub, *lrl = NULL, *lru = NULL;
    int ii, *rowInd, no = 0;
    mwIndex startRow, stopRow; 
    char msgbuf[1024];
    
    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    mwIndex *H_ir, *H_jc;
    int *rows = NULL;
    CoinBigIndex *cols = NULL;
    
    //SOS
    CbcObject ** objects;
    int no_sets;
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
	f = mxGetPr(prhs[0]);
	A = mxGetPr(prhs[1]); 
    A_ir = mxGetIr(prhs[1]);
    A_jc = mxGetJc(prhs[1]);
    rl = mxGetPr(prhs[2]);
    ru = mxGetPr(prhs[3]);
    lb = mxGetPr(prhs[4]); 
    ub = mxGetPr(prhs[5]);
    if(mxIsInt32(prhs[6])) //backwards compatibility
        ixtype = (int*)mxGetData(prhs[6]);
    else
        xtype = mxArrayToString(prhs[6]);
    if(nrhs > 9) //optional quadratic part
        H = mxGetPr(prhs[9]);

    //Get sizes
    ndec = mxGetM(prhs[0]);
    ncon = mxGetM(prhs[1]);   
    
    //Get Options if Specified
    if(nrhs > 8) {
    	if(mxGetField(prhs[8],0,"intTol"))
            intTol = *mxGetPr(mxGetField(prhs[8],0,"intTol"));
        if(mxGetField(prhs[8],0,"maxnodes"))
            maxnodes = (int)*mxGetPr(mxGetField(prhs[8],0,"maxnodes"));
        if(mxGetField(prhs[8],0,"maxtime"))
            maxtime = (int)*mxGetPr(mxGetField(prhs[8],0,"maxtime"));
        if(mxGetField(prhs[8],0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(prhs[8],0,"display"));
        if(mxGetField(prhs[8],0,"debug"))
            debug = (int)*mxGetPr(mxGetField(prhs[8],0,"debug"));
        if(mxGetField(prhs[8],0,"objbias"))
            objc = *mxGetPr(mxGetField(prhs[8],0,"objbias"));
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
    contobj = mxGetPr(plhs[4]);
    
    try
    {
        //Objects
        CbcModel cbcmodel;
        CoinModel model;
        OsiClpSolverInterface OSImodel;
        DerivedHandler *mexprinter;     
        DerivedEvent *ctrlCEvent;
        
        //Allocate Index Vector
        rowInd = (int*)mxCalloc(ncon,sizeof(int)); //set as max no of rows
        
        //Process Bounds (must copy in as we will change them)
        llb = (double*)mxCalloc(ndec,sizeof(double));
        lub = (double*)mxCalloc(ndec,sizeof(double));   
        //Create bounds if empty
        if(mxIsEmpty(prhs[4])) {
            for(i=0;i<ndec;i++)
                llb[i] = -COIN_DBL_MAX;
        }
        else
            memcpy(llb,lb,ndec*sizeof(double));

        if(mxIsEmpty(prhs[5])) {
            for(i=0;i<ndec;i++)
                lub[i] = COIN_DBL_MAX;
        }
        else
            memcpy(lub,ub,ndec*sizeof(double));

        //Ensure 'finite' bounds
        for(i = 0; i < ndec; i++) {
            if(mxIsInf(llb[i]))
                llb[i] = -COIN_DBL_MAX;
            if(mxIsInf(lub[i]))
                lub[i] = COIN_DBL_MAX;        
        }
        
        //Modify bounds based on binary constraints
        if(mxIsChar(prhs[6])) {
            for(ii = 0; ii < ndec; ii++) {
                switch(xtype[ii])
                {
                    case 'B':
                        //Enforce binary bounds if not specified
                        if(llb[ii] == -COIN_DBL_MAX)
                            llb[ii] = 0;
                        if(lub[ii] == COIN_DBL_MAX)
                            lub[ii] = 1;
                        break;
                }
            }
        }       

        //Add Linear Constraints (Sparse)
        if(ncon) {
            for(i = 0; i < ndec; i++) {
                startRow = A_jc[i];
                stopRow = A_jc[i+1];
                no = (int)(stopRow - startRow);
                if(no > 0) {
                    for(j = 0, k = startRow; k < stopRow; j++, k++) //build int32 row indicies
                        rowInd[j] = (int)A_ir[k];                
                }
                model.addColumn(no,rowInd,&A[startRow],llb[i],lub[i],f[i]);
            }
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
        else {//just bounds
            for(ii=0;ii<ndec;ii++) {
                model.setObjective(ii,f[ii]);
                model.setColumnBounds(ii,llb[ii],lub[ii]);
            }
        }
        
        //Add Row Bounds
        for (ii = 0; ii < ncon; ii++)
        	model.setRowBounds(ii, lrl[ii], lru[ii]);
          
        //Add Objective Offset
        model.setObjectiveOffset(-objc);
        
        //Add Integer and Binary Constraints
        for(ii = 0; ii < ndec; ii++) {
            if(mxIsInt32(prhs[6])) {
                if(ixtype[ii])
                    model.setInteger(ii);
            }
            else {
                switch(xtype[ii])
                {
                    case 'C': break;
                    case 'I': 
                    case 'B':
                        model.setInteger(ii); 
                        break;
                    default:
                        throw exception("Unknown xtype, only 'C', 'I' and 'B' are accepted");
                }
            }
        } 
        
        //Add Quadratic Objective (can't work this out - suggestions anyone??)
        if(H && !mxIsEmpty(prhs[9])) {          
//            H_ir = mxGetIr(prhs[7]);
//            H_jc = mxGetJc(prhs[7]);
//            //Convert Indicies
//            mwIndex nzH = H_jc[ndec];
//            rows = (int*)mxCalloc(nzH,sizeof(int));
//            cols = (CoinBigIndex*)mxCalloc(ndec+1,sizeof(CoinBigIndex));
//            //Assign Convert Data Type Vectors
//            for(i = 0; i <= ndec; i++)
//                cols[i] = (CoinBigIndex)H_jc[i];
//            for(i = 0; i < nzH; i++)
//                rows[i] = (int)H_ir[i];
//            Load QuadObj (don't know what to do here...)
//            clpQuad = new ClpQuadraticObjective(f, ndec, cols, rows, H);
        }
        
        //Load Problem into OSI Interface
        OSImodel.loadFromCoinModel(model);
        OSImodel.setHintParam(OsiDoReducePrint,true,OsiHintTry); //prevent clp output
        //Load into CBC
        cbcmodel = CbcModel(OSImodel);        
                
        //Load Cutting Generators + Heuristics (Following examples .. not really sure here!)
        CglProbing generator1;
        generator1.setUsingObjective(true);
        generator1.setMaxPass(1);
        generator1.setMaxPassRoot(5);
        // Number of unsatisfied variables to look at
        generator1.setMaxProbe(10);
        generator1.setMaxProbeRoot(1000);
        // How far to follow the consequences
        generator1.setMaxLook(50);
        generator1.setMaxLookRoot(500);
        // Only look at rows with fewer than this number of elements
        generator1.setMaxElements(200);
        generator1.setRowCuts(3);
        //Other Generators
        CglMixedIntegerRounding2 mixedGen;
        CglFlowCover flowGen;        
        // Add in generators
        cbcmodel.addCutGenerator(&generator1,-1,"Probing");
        cbcmodel.addCutGenerator(&flowGen,-1,"FlowCover");
        cbcmodel.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
        //Heuristics
        CbcRounding heuristic1(cbcmodel);
        cbcmodel.addHeuristic(&heuristic1);
        CbcHeuristicLocal heuristic2(cbcmodel);
        cbcmodel.addHeuristic(&heuristic2);
        
        //Setup SOS
        if(nrhs > 7 && !mxIsEmpty(prhs[7])) {
            no_sets = (int)mxGetNumberOfElements(mxGetField(prhs[7],0,"type"));
            if(no_sets > 0) {
                //Collect Types
                sostype = mxArrayToString(mxGetField(prhs[7],0,"type"));
                //Allocate Set Memory
                objects = new CbcObject * [no_sets];
                //Copy in SOS data, creating CbcSOS objects as we go
                for(i=0;i<no_sets;i++) {
                    type = (int)sostype[i] - '0'; //convert to numerical representation
                    if(type < 1 || type > 2)
                        throw exception("Only SOS of type '1' and '2' are supported");
                    no_entries = (int)mxGetNumberOfElements(mxGetCell(mxGetField(prhs[7],0,"index"),i));
                    //Create sosind memory and copy in
                    isosind = new int[no_entries];
                    sosind = mxGetPr(mxGetCell(mxGetField(prhs[7],0,"index"),i));
                    for(j=0;j<no_entries;j++)
                        isosind[j] = (int)sosind[j]-1;                
                    //Get soswt
                    soswt = mxGetPr(mxGetCell(mxGetField(prhs[7],0,"weight"),i));
                    //Create Set object
                    objects[i] = new CbcSOS(&cbcmodel,no_entries,isosind,soswt,(int)i,type);
                    delete isosind; //free memory as we go round, CoinSet copies internally
                }

                //Add objects to model
                cbcmodel.addObjects(no_sets,objects); 
                //Priorities??

                //Delete objects
                for(i=0;i<no_sets;i++)
                    delete objects[i];
                delete [] objects;
                //Delete type string
                mxFree(sostype);
            }
        }
        
        //Set Options
        cbcmodel.setMaximumNodes(maxnodes);
		cbcmodel.setMaximumSeconds(maxtime);
        cbcmodel.setIntegerTolerance(intTol);
        if(printLevel) {
            mexprinter = new DerivedHandler();
            mexprinter->setLogLevel(0,printLevel);
            cbcmodel.passInMessageHandler(mexprinter);
        }
        //Add Event Handler for Ctrl+C
        ctrlCEvent = new DerivedEvent();      
        cbcmodel.passInEventHandler(ctrlCEvent); 
        
        //Solve using CBC + CLP
        cbcmodel.initialSolve();    //relaxed LP
        cbcmodel.branchAndBound();  //full solve
        
        //Assign Return Arguments
        sol = cbcmodel.getColSolution();
        
        if(sol != NULL) {
            memcpy(x,sol,ndec*sizeof(double));
            *fval = cbcmodel.getObjValue();
            *exitflag = getStatus(cbcmodel.secondaryStatus());
            *iter = cbcmodel.getNodeCount();
            *contobj = cbcmodel.getContinuousObjective();
        }

        //Clean up memory
        mxFree(rowInd);
        mxFree(llb);
        mxFree(lub);
        if(xtype != NULL) mxFree(xtype); xtype = NULL;
        if(lrl != NULL) mxFree(lrl); lrl = NULL;
        if(lru != NULL) mxFree(lru); lru = NULL;
        if(printLevel)
            delete mexprinter;
        if(rows) mxFree(rows);
        if(cols) mxFree(cols);
        
    }
    //Error Handling (still crashes Matlab though...)
    catch(CoinError e)
    {
        sprintf(msgbuf,"Caught Coin Error: %s",e.message());
        mexErrMsgTxt(msgbuf);
    }
    catch(exception& e)
    {
        sprintf(msgbuf,"Caught CBC Error: %s",e.what());
        mexErrMsgTxt(msgbuf);           
    }  
}               


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, ncon;
    
    //Correct number of inputs
    if(nrhs < 7)
        mexErrMsgTxt("You must supply at least 7 arguments to cbc (f, A, rl, ru, lb, ub, xtype)"); 
    
    //Check we have an objective
    if(mxIsEmpty(prhs[0]))
        mexErrMsgTxt("You must supply an objective function!");
    
    //Check we have some constraints
    if(mxIsEmpty(prhs[1]) && mxIsEmpty(prhs[4]) && mxIsEmpty(prhs[5]))
        mexErrMsgTxt("You have not supplied any constraints!");
   
    //Check SOS structure
    if(nrhs > 7 && !mxIsEmpty(prhs[7])) {
        if(!mxIsStruct(prhs[7]))
            mexErrMsgTxt("The SOS argument must be a structure!");       
        if(mxGetFieldNumber(prhs[7],"type") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'type'");
        if(mxGetFieldNumber(prhs[7],"index") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'index'");
        if(mxGetFieldNumber(prhs[7],"weight") < 0)
            mexErrMsgTxt("The sos structure should contain the field 'weight'");
        //Ensure type is char array
        if(!mxIsChar(mxGetField(prhs[7],0,"type")))
            mexErrMsgTxt("sos.type should be a char array");
        //Check multiple sets length
        int no_sets = (int)mxGetNumberOfElements(mxGetField(prhs[7],0,"type")); 
        if(no_sets > 1) {
            if(!mxIsCell(mxGetField(prhs[7],0,"index")) || mxIsEmpty(mxGetField(prhs[7],0,"index")))
                mexErrMsgTxt("sos.index must be a cell array, and not empty!");
            if(!mxIsCell(mxGetField(prhs[7],0,"weight")) || mxIsEmpty(mxGetField(prhs[7],0,"weight")))
                mexErrMsgTxt("sos.weight must be a cell array, and not empty!");
            if(mxGetNumberOfElements(mxGetField(prhs[7],0,"index")) != no_sets)
                mexErrMsgTxt("sos.index cell array is not the same length as sos.type!");
            if(mxGetNumberOfElements(mxGetField(prhs[7],0,"weight")) != no_sets)
                mexErrMsgTxt("sos.weight cell array is not the same length as sos.type!");        
        }
    }
    
    //Check options is a structure with correct fields
    if(nrhs > 8 && !mxIsStruct(prhs[8]))
        mexErrMsgTxt("The options argument must be a structure!");
    
    //Get Sizes
    ndec = mxGetM(prhs[0]);
    ncon = mxGetM(prhs[1]);
    
    //Check xtype
    if(mxIsEmpty(prhs[6]) || mxGetNumberOfElements(prhs[6]) < ndec || (!mxIsChar(prhs[6]) && !mxIsInt32(prhs[6])))
        mexErrMsgTxt("xtype should be a char array with ndec elements");
    
    //Check Constraint Pairs
    if(ncon && mxIsEmpty(prhs[2]))
        mexErrMsgTxt("rl is empty!");
    if(ncon && mxIsEmpty(prhs[3]))
        mexErrMsgTxt("ru is empty!");
    
    //Check Sparsity (only supported in A and H)
    if(mxIsSparse(prhs[0]))
        mexErrMsgTxt("Only A is a sparse matrix");
    if(!mxIsSparse(prhs[1]))
        mexErrMsgTxt("A must be a sparse matrix");
    if(nrhs > 9 && !mxIsSparse(prhs[9]))
        mexErrMsgTxt("H must be a sparse matrix");
    
     //Check Orientation
    if(mxGetM(prhs[0]) < mxGetN(prhs[0]))
        mexErrMsgTxt("f must be a column vector");
    if(mxGetM(prhs[2]) < mxGetN(prhs[2]))
        mexErrMsgTxt("rl must be a column vector");
    if(mxGetM(prhs[3]) < mxGetN(prhs[3]))
        mexErrMsgTxt("ru must be a column vector");
    if(mxGetM(prhs[4]) < mxGetN(prhs[4]))
        mexErrMsgTxt("lb must be a column vector");
    if(mxGetM(prhs[5]) < mxGetN(prhs[5]))
        mexErrMsgTxt("ub must be a column vector");
    
    //Check Sizes
    if(ncon) {
        if(mxGetN(prhs[1]) != ndec)
            mexErrMsgTxt("A has incompatible dimensions");
        if(mxGetM(prhs[2]) != ncon)
            mexErrMsgTxt("rl has incompatible dimensions");
        if(mxGetM(prhs[3]) != ncon)
            mexErrMsgTxt("ru has incompatible dimensions");
    }
    if(!mxIsEmpty(prhs[4]) && (mxGetM(prhs[4]) != ndec))
        mexErrMsgTxt("lb has incompatible dimensions");
    if(!mxIsEmpty(prhs[5]) && (mxGetM(prhs[5]) != ndec))
        mexErrMsgTxt("ub has incompatible dimensions");    
    if(nrhs > 9 && !mxIsEmpty(prhs[9]) && ((mxGetM(prhs[9]) != ndec) || (mxGetN(prhs[9]) != ndec)))
        mexErrMsgTxt("H has incompatible dimensions");
}


//Convert exiflag to OPTI status
double getStatus(int stat)
{
    double ret = -1;

    switch(stat)
    {
        case 0: //looks optimal
            ret = 1;
            break;
        case 1: //lp relax infeasible
            ret = -1;
            break;
        case 2: //gap reached
        case 3: //max nodes
        case 4: //max time
            ret = 0;
            break;
        case 5: //user exit
            ret = -5;
            break;
        default: //inaccuracy / unbounded
            ret = -2;
            break;
    }
    return ret;
}  

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" CBC: COIN-OR Branch and Cut [v%s]\n",CBC_VERSION);
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/Cbc\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - CGL    [v%s] (Eclipse Public License)\n",CGL_VERSION);
    mexPrintf("  - CLP    [v%s] (Eclipse Public License)\n",CLP_VERSION);
    mexPrintf("  - CoinUtils [v%s] (Eclipse Public License)\n",COINUTILS_VERSION);
    mexPrintf("  - OSI    [v%s] (Eclipse Public License)\n",OSI_VERSION);
    
    mexPrintf("\n MEX Interface J.Currie 2013 (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}