/* SCIPSDPMEX - A MATLAB MEX Interface to SCIP-SDP
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2014
 * www.i2c2.aut.ac.nz
 */

#include "mex.h"
#include <exception>
#include <ctype.h>
#include <stdio.h>

#include "SdpCone.h"
#include "objconshdlr_sdp.h"
#include "objrelax_sdp.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "spxdefines.h"
#include "scipmex.h"
#include "mkl.h"
#include "dsdp5.h"

//Enable for Debug print out
#define DEBUG

//DSDP Version
#define DSDP_VERSION "5.8"
#define SCIPSDP_VERSION "1.0"

using namespace std;

//Argument Enumeration (in expected order of arguments)
enum {eF, eA, eB, eLB, eUB, eSDP, eXTYPE, eOPTS};
//PRHS Defines    
#define pF      prhs[eF]
#define pA      prhs[eA]
#define pB      prhs[eB]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pSDP    prhs[eSDP]
#define pXTYPE  prhs[eXTYPE]
#define pOPTS   prhs[eOPTS]

//Function Prototypes
void printSolverInfo();
void addSDPCone(SCIP* scip, SCIP_VAR **scipvars, const mxArray *cone, int block);
void checkInputs(const mxArray *prhs[], int nrhs);
void getIntOption(const mxArray *opts, const char *option, int &var);
void getDblOption(const mxArray *opts, const char *option, double &var);
void getStrOption(const mxArray *opts, const char *option, char *str);

//Message Handler Callback
void msginfo(SCIP_MESSAGEHDLR *messagehdlr, FILE *file, const char *msg)
{
    mexPrintf(msg);
    mexEvalString("drawnow;"); //flush draw buffer
}
//Message Buffer
char msgbuf[BUFSIZE];

//Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *f, *A, *b, *lb, *ub;
    char *xtype;
    
    //Return Args
    double *x, *fval, *exitflag, *iter, *nodes, *gap;
    const char *fnames[3] = {"LPiter","BBnodes","BBgap"};
    
    //Common Options
    int maxiter = 1500;
    int maxnodes = 10000;
    int maxpresolve = -1;
    double maxtime = 1000;
    double primtol = 1e-6;
    double objbias = 0.0;
    int printLevel = 0;
    
    //Internal Vars
    size_t ncones = 0, ndec = 0, ncon = 0;
    size_t ncnt = 0, nint = 0, nbin = 0;
    size_t i, j; 
    int alb = 0, aub = 0, no = 0;
        
    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    mwIndex startRow, stopRow;
    
    //Print Header or return version string if requested
    if(nrhs < 1) {
        if(nlhs < 1) 
            printSolverInfo();
        else {
            plhs[0] = mxCreateString(SCIPSDP_VERSION);
        }
        return;
    }        
    
    //Check Inputs
    checkInputs(prhs,nrhs); 

    //Get pointers to input vars
    f = mxGetPr(pF);
    A = mxGetPr(pA); A_ir = mxGetIr(pA); A_jc = mxGetJc(pA);
    b = mxGetPr(pB);
    lb = mxGetPr(pLB); ub = mxGetPr(pUB);
    if(nrhs > eSDP && !mxIsEmpty(pSDP)) {
        if(mxIsCell(pSDP))
            ncones = mxGetNumberOfElements(pSDP);
        else
            ncones = 1;
    }
    if(nrhs > eXTYPE)
        xtype = mxArrayToString(pXTYPE);
    //Get sizes from input args
    ndec = mxGetNumberOfElements(pF);
    ncon = mxGetM(pA); 
                        
    //Get Common Options if specified
    if(nrhs > eOPTS) {
        getIntOption(pOPTS,"maxiter",maxiter);
        getIntOption(pOPTS,"maxnodes",maxnodes);
        getIntOption(pOPTS,"maxpresolve",maxpresolve);
        getDblOption(pOPTS,"maxtime",maxtime);
        getDblOption(pOPTS,"tolrfun",primtol);
        getDblOption(pOPTS,"objbias",objbias);
        getIntOption(pOPTS,"display",printLevel);
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
    //Statistic Structure Output
    plhs[3] = mxCreateStructMatrix(1,1,3,fnames);
    mxSetField(plhs[3],0,fnames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[2],mxCreateDoubleMatrix(1,1, mxREAL));
    iter  = mxGetPr(mxGetField(plhs[3],0,fnames[0]));
    nodes = mxGetPr(mxGetField(plhs[3],0,fnames[1]));
    gap   = mxGetPr(mxGetField(plhs[3],0,fnames[2]));

    //SCIP Objects
    SCIP* scip;
    SCIP_VAR** vars = NULL;
    SCIP_CONS** cons = NULL;
    SCIP_VAR *objb = NULL;

    //Create SCIP Object
    SCIP_ERR( SCIPcreate(&scip) , "Error creating SCIP object");
    //Add SCIP-SDP Plugins
    SCIP_ERR( SCIPincludeObjConshdlr(scip, new ObjConshdlrSdp(scip), TRUE), "Error including SCIP-SDP Constraint Handler plugin");
    SCIP_ERR( SCIPincludeObjRelax(scip, new scip::ObjRelaxSdp(scip), TRUE), "Error including SCIP-SDP Constraint Handler plugin");
    const char *name = "sdpsolver";
    const char *desc = "which sdpsolver should be called";
    SCIP_PARAMDATA *paramdata = NULL;
    SCIP_ERR( SCIPaddStringParam(scip, name, desc, NULL, FALSE, "dsdp" , NULL, paramdata), "Error adding SDP Solver Parameter");	
    //Add default plugins
    SCIP_ERR( SCIPincludeDefaultPlugins(scip), "Error including SCIP default plugins");
    //Add Ctrl-C Event Handler
    SCIP_ERR( SCIPincludeCtrlCEventHdlr(scip), "Error adding Ctrl-C Event Handler");
    
    //Create Empty Problem
    SCIP_ERR( SCIPcreateProbBasic(scip,"OPTI Problem"), "Error creating basic SCIP problem");    

    //Create continuous xtype array if empty or not supplied
    if(nrhs <= eXTYPE || mxIsEmpty(pXTYPE)) {
        xtype = (char*)mxCalloc(ndec,sizeof(char));
        for(i=0;i<ndec;i++)
            xtype[i] = 'c';
    }

    //Create infinite bounds if empty
    if(mxIsEmpty(pLB)) {
        lb = (double*)mxCalloc(ndec,sizeof(double)); alb=1;
        for(i=0;i<ndec;i++)
            lb[i] = -1e50;
    }
    if(mxIsEmpty(pUB)) {
        ub = (double*)mxCalloc(ndec,sizeof(double)); aub=1;
        for(i=0;i<ndec;i++)
            ub[i] = 1e50;
    }

    //Create SCIP Variables (also loads linear objective + bounds)
    SCIP_ERR( SCIPallocMemoryArray(scip,&vars,(int)ndec), "Error allocating variable memory");
    double llb, lub;
    for(i=0;i<ndec;i++)
    {
        SCIP_VARTYPE vartype;
        //Assign variable type
        switch(tolower(xtype[i]))
        {
            case 'i':
                vartype = SCIP_VARTYPE_INTEGER; 
                llb = lb[i]; lub = ub[i]; 
                sprintf(msgbuf,"ivar%d",nint++);
                break;
            case 'b':
                vartype = SCIP_VARTYPE_BINARY; 
                llb = lb[i] <= -1e50 ? 0 : lb[i]; //if we don't do this, SCIP fails during presolve
                lub = ub[i] >= 1e50 ? 1 : ub[i];
                sprintf(msgbuf,"bvar%d",nbin++);
                break;
            case 'c':
                vartype = SCIP_VARTYPE_CONTINUOUS; 
                llb = lb[i]; lub = ub[i]; 
                sprintf(msgbuf,"xvar%d",ncnt++);
                break;
            default:
                sprintf(msgbuf,"Unknown variable type for variable %d",i);
                mexErrMsgTxt(msgbuf);
        }
        //Create variable
        SCIP_ERR( SCIPcreateVarBasic(scip,&vars[i],msgbuf,llb,lub,f[i],vartype), "Error creating basic SCIP variable");
        //Add to problem
        SCIP_ERR( SCIPaddVar(scip,vars[i]), "Error adding SCIP variable to problem");
    }

    //Add objective bias term if non-zero
    if(objbias != 0) {
        SCIP_ERR( SCIPcreateVarBasic(scip, &objb, "objbiasterm", objbias, objbias, 1.0, SCIP_VARTYPE_CONTINUOUS), "Error adding objective bias variable");
        SCIP_ERR( SCIPaddVar(scip, objb), "Error adding objective bias variable");
    }

    //Add Linear Constraints (if they exist)
    if(ncon) {     
        //Allocate memory for all constraints (we create them all now, as we have to add coefficients in column order)
        SCIP_ERR( SCIPallocMemoryArray(scip, &cons, (int)ncon), "Error allocating constraint memory");
        //Create each constraint and add row bounds, but leave coefficients empty
        for(i=0;i<ncon;i++) {
            SCIPsnprintf(msgbuf, BUFSIZE, "lincon%d", i); //appears constraints require a name
            SCIP_ERR( SCIPcreateConsBasicLinear(scip,&cons[i],msgbuf,0,NULL,NULL,-SCIPinfinity(scip),b[i]), "Error creating basic SCIP linear constraint");
        }
        //Now for each column (variable), add coefficients
        for(i = 0; i < ndec; i++) {
            //Determine number of nz in this column
            startRow = A_jc[i];
            stopRow = A_jc[i+1];
            no = (int)(stopRow - startRow);
            //If we have nz in this column
            if(no > 0) {
                //Add each coefficient
                for(j = startRow; j < stopRow; j++)
                    SCIP_ERR( SCIPaddCoefLinear(scip, cons[A_ir[j]], vars[i], A[j]), "Error adding constraint linear coefficient");
            }
        }
        //Now for each constraint, add it to the problem, then release it
        for(i=0;i<ncon;i++) {
            SCIP_ERR( SCIPaddCons(scip,cons[i]), "Error adding linear constraint");                   
            SCIP_ERR( SCIPreleaseCons(scip,&cons[i]), "Error releasing linear constraint");
        }
    }
 
    //Add Semidefinite Constraints
    for(i=0;i<ncones;i++) {
        if(ncones == 1 && !mxIsCell(pSDP)) 
            addSDPCone(scip,vars,pSDP,(int)i);
        else 
            addSDPCone(scip,vars,mxGetCell(pSDP,i),(int)i);
    }
    
    //Set SCIP-SDP Recommended Options
    SCIP_ERR( SCIPsetIntParam(scip,"relaxing/SDPRelax/freq",0), "Error setting SDPrelaxfreq");
    SCIP_ERR( SCIPsetIntParam(scip, "lp/solvefreq", 1), "Error setting lp/solvefreq");
    SCIP_ERR( SCIPsetRealParam(scip, "numerics/epsilon", 1e-6), "Error setting numerics/epsilon" );
    SCIP_ERR( SCIPsetRealParam(scip, "numerics/feastol", 1e-4), "Error setting numerics/feastol");
    SCIP_ERR( SCIPsetStringParam(scip, "sdpsolver", "dsdp"), "Error setting sdpsolver");
    SCIP_ERR( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE), "Error setting lp/cleanuprows");
    SCIP_ERR( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE), "Error setting lp/cleanuprowsroot");
    SCIP_ERR( SCIPsetIntParam(scip, "lp/rowagelimit", 10), "Error setting lp/rowagelimit");
    SCIP_ERR( SCIPsetIntParam(scip, "separating/cutagelimit", 10), "Error setting separating/cutagelimit");
    SCIP_ERR( SCIPsetIntParam(scip, "separating/maxrounds", 20), "Error setting separating/maxrounds");
    SCIP_ERR( SCIPsetIntParam(scip, "separating/intobj/freq", -1), "Error setting separating/intobj/freq");
    SCIP_ERR( SCIPsetIntParam(scip, "branching/inference/priority", -500000000), "Error setting branching/inference/priority"); //turn off basically

    //Set Common OPTI Options  
    SCIP_ERR( SCIPsetRealParam(scip,"limits/time",maxtime), "Error setting maxtime");
    SCIP_ERR( SCIPsetLongintParam(scip,"lp/iterlim",maxiter), "Error setting iterlim");
    SCIP_ERR( SCIPsetLongintParam(scip,"limits/nodes",maxnodes), "Error setting nodes");
    SCIP_ERR( SCIPsetRealParam(scip,"numerics/lpfeastol",primtol), "Error setting lpfeastol"); 
    SCIP_ERR( SCIPsetIntParam(scip,"presolving/maxrounds",maxpresolve), "Error setting max presolve rounds"); 
    
    //If user has requested print out
    if(printLevel)
    {
        //Create Message Handler
        SCIP_MESSAGEHDLR *mexprinter;
        SCIPmessagehdlrCreate(&mexprinter,TRUE,NULL,FALSE,&msginfo,&msginfo,&msginfo,NULL,NULL);
        SCIP_ERR( SCIPsetMessagehdlr(scip,mexprinter), "Error adding message handler");
        //Set Verbosity Level
        SCIP_ERR( SCIPsetIntParam(scip,"display/verblevel",printLevel), "Error setting verblevel");
    }
    
    //Solve Problem
    SCIP_ERR( SCIPsolve(scip), "Error solving SCIP problem!");
    //Assign Return Arguments
    if(SCIPgetNSols(scip) > 0 ) 
    {
        SCIP_SOL* scipbestsol = SCIPgetBestSol(scip);
        //Assign x
        for(i = 0;i<ndec;i++)
           x[i] = SCIPgetSolVal(scip,scipbestsol,vars[i]);            
        //Assign fval
        *fval = SCIPgetSolOrigObj(scip, scipbestsol);
        //Get Solve Statistics
        *iter = (double)SCIPgetNLPIterations(scip);
        *nodes = (double)SCIPgetNTotalNodes(scip);
        *gap = SCIPgetGap(scip);
    }
    //Get Solution Status
    *exitflag = (double)SCIPgetStatus(scip);

    //Clean up memory from MATLAB mode*/
    mxFree(xtype);
    if(alb) mxFree(lb); alb = 0;
    if(aub) mxFree(ub); aub = 0;

    //Release Variables
    for(i=0;i<ndec;i++)
        SCIP_ERR( SCIPreleaseVar(scip,&vars[i]), "Error releasing SCIP variable");
    if(objb != NULL)
        SCIP_ERR( SCIPreleaseVar(scip,&objb), "Error releasing SCIP objective bias variable");    
    
    //Now free SCIP arrays & problem
    SCIPfreeMemoryArray(scip, &vars);
    if(ncon) SCIPfreeMemoryArray(scip, &cons);
          
    //Clean up general SCIP memory
    SCIP_ERR( SCIPfree(&scip), "Error releasing SCIP problem");
}      

//Add SDP Constraint
void addSDPCone(SCIP* scip, SCIP_VAR **scipvars, const mxArray *cone, int block)
{
    size_t i,j,idx,midx;
    double *SDP_pr  = mxGetPr(cone);
    mwIndex *SDP_ir = mxGetIr(cone);    
    mwIndex *SDP_jc = mxGetJc(cone);
    int SDP_M       = (int)mxGetM(cone);
    int SDP_N       = (int)mxGetN(cone); //remember [C A0 A1 A2...] so non-square    
    int SDP_C_nnz   = 0; //nnz in C
    int SDP_A_nnz   = 0; //nnz in current A 
    int SDP_DIM     = (int)(sqrt((double)SDP_M)); //calculate dimension
    int rind, cind;
    
    //Find NNZ
    SDP_C_nnz = (int)(SDP_jc[1]-SDP_jc[0]);
    SDP_A_nnz = (int)SDP_jc[SDP_N] - SDP_C_nnz;
    
    #ifdef DEBUG
        mexPrintf("SDP_DIM [block %d]: %d, M: %d, N: %d\n",block,SDP_DIM,SDP_M,SDP_N);
        mexPrintf("C nnz: %d, ALL A nnz: %d\n",SDP_C_nnz,SDP_A_nnz);
    #endif
        
    //Allocate constraint memory (we allocate too much here...)
    SCIP_VAR ** vars = NULL;
    int *col = NULL, *row = NULL, *const_col = NULL, *const_row = NULL, nnza = 0, nnzc = 0;
    double *vals = NULL, *const_vals = NULL;            
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &vars, SDP_A_nnz), "Error Allocating SCIP-SDP Variable Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &col, SDP_A_nnz), "Error Allocating SCIP-SDP Column Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &row, SDP_A_nnz), "Error Allocating SCIP-SDP Row Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &vals, SDP_A_nnz), "Error Allocating SCIP-SDP Values Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &const_col, SDP_C_nnz), "Error Allocating SCIP-SDP Constant Column Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &const_row, SDP_C_nnz), "Error Allocating SCIP-SDP Constant Row Memory");
    SCIP_ERR(SCIPallocBlockMemoryArray(scip, &const_vals, SDP_C_nnz), "Error Allocating SCIP-SDP Constant Value Memory");    
        
    //Copy in C
    idx = 0; midx = 0;
    for(i=0;i<SDP_C_nnz;i++) {
        //Row & Col Index
        rind = SDP_ir[i] % SDP_DIM;
        cind = (int)((SDP_ir[i] - rind)/SDP_DIM);
        if(rind <= cind) {
            //Copy In
            const_row[idx] = rind+1;
            const_col[idx] = cind+1;
            const_vals[idx] = SDP_pr[i];
            nnzc++;            
            #ifdef DEBUG
                mexPrintf("(%d) - C[%d,%d] = %f\n",idx,const_row[idx],const_col[idx],const_vals[idx]);
            #endif
            idx++;
        }
        midx++;
    } 
    idx = 0;
    //Copy in all As
    for(i=1;i<SDP_N;i++) { //for each variable (A_i matrix)
        for(j=SDP_jc[i];j<SDP_jc[i+1];j++) { //for each element in A_i
            //Row & Col Index
            rind = SDP_ir[midx] % SDP_DIM;
            cind = (int)((SDP_ir[midx] - rind)/SDP_DIM);
            if(rind <= cind) {
                //Copy in
                vars[idx] = scipvars[i-1];
                row[idx] = rind+1;
                col[idx] = cind+1;
                vals[idx] = SDP_pr[midx];
                nnza++;                
                #ifdef DEBUG
                    mexPrintf("(%d) - A[%d][%d,%d] = %f\n",idx,i,row[idx],col[idx],vals[idx]);
                #endif   
                idx++;
            }            
            midx++;
        }         
    }
    
    //Create SCIP-SDP Constraint
    SdpCone sdpcone(scip, SDP_DIM, vars, col, row, vals, nnza, const_col, const_row, const_vals, nnzc);
    #ifdef DEBUG
        mexPrintf("Actual NNZC %d, NNZA %d\n",nnzc,nnza);
        mexPrintf("Created Sdcone %d\n",block);
    #endif 
    //Create SCIP Constraint
    SCIP_CONS* sdpcon;
    SCIPsnprintf(msgbuf, BUFSIZE, "SDP-Constraint-%d", block);
    SCIP_ERR( SCIPcreateConsSdp(scip, &sdpcon, msgbuf, sdpcone), "Error Creating SDP Constraint" );
    SCIP_ERR( SCIPaddCons(scip, sdpcon), "Error Adding SDP Constraint" );
    SCIP_ERR( SCIPreleaseCons(scip, &sdpcon), "Error Releasing SDP Constraint" );
    #ifdef DEBUG
        mexPrintf("Added Sdcone %d\n",block);
    #endif 
    /*//Copy cast ir
    iSDP_ir = (int*)mxCalloc(SDP_jc[SDP_N],sizeof(int));
    for(i=0;i<(int)SDP_jc[SDP_N];i++)
        iSDP_ir[i] = (int)SDP_ir[i];

    //Set Cone Size
    DSDP_ERR( SDPConeSetBlockSize(sdpcone,block,(int)SDP_DIM), "Error setting cone dimensions"); 
    DSDP_ERR( SDPConeUsePackedFormat(sdpcone, block), "Error setting cone packed format");
    DSDP_ERR( SDPConeSetSparsity(sdpcone,block,(int)SDP_jc[SDP_N]), "Error setting cone sparsity nnz");
    //Add C
    DSDP_ERR( SDPConeSetASparseVecMat(sdpcone,block,0,SDP_DIM,1.0,0,iSDP_ir,SDP_pr,SDP_C_nnz), "Error setting cone C matrix");
    //Add Each A
    SDP_A_index = SDP_C_nnz; //set index to elements after C
    for(i=1;i<SDP_N;i++) {
        SDP_A_nnz = (int)(SDP_jc[i+1]-SDP_jc[i]);  
        #if DEBUG
            mexPrintf("A[%d] nnz: %d, A index: %d, pr[0] =  %f ir[0] = %d\n",i-1,SDP_A_nnz,SDP_A_index,SDP_pr[SDP_A_index],iSDP_ir[SDP_A_index]);
        #endif
        DSDP_ERR( SDPConeSetASparseVecMat(sdpcone,block,i,SDP_DIM,1.0,0,&iSDP_ir[SDP_A_index],&SDP_pr[SDP_A_index],SDP_A_nnz), "Error setting cone A matrix");
        SDP_A_index += SDP_A_nnz; //shift index forward
    }
    //Do not free memory here, MATLAB will take care of it once MEX completes (DSDP will use it during solving)  [See MATLAB/User's Guide/C/C++ and Fortran API Reference/mxFree] */
}


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, ncon;    
    
    //Correct number of inputs
    if(nrhs <= eUB)
        mexErrMsgTxt("You must supply at least 6 arguments to scipsdp (f, A, b, lb, ub, sdcone)"); 
    
    //Check we have an objective
    if(mxIsEmpty(pF))
        mexErrMsgTxt("You must supply a linear objective function via f (all zeros if not required)!");
    
    //Check options is a structure
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && !mxIsStruct(pOPTS))
        mexErrMsgTxt("The options argument must be a structure!");
    
    //Get Sizes    
    ndec = mxGetNumberOfElements(pF);
    ncon = mxGetM(pA);
    
    //Check Constraint Pairs
    if(ncon && mxIsEmpty(pB))
        mexErrMsgTxt("b is empty!");

    //Check Sparsity (only supported in A and H)
    if(!mxIsEmpty(pA)) {
        if(mxIsSparse(pF) || mxIsSparse(pB) || mxIsSparse(pLB))
            mexErrMsgTxt("Only A is a sparse matrix");
        if(!mxIsSparse(pA))
            mexErrMsgTxt("A must be a sparse matrix");
    }
    
    //Check xtype data type
    if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE) && mxGetClassID(pXTYPE) != mxCHAR_CLASS)
        mexErrMsgTxt("xtype must be a char array");
  
    
    //Check Sizes
    if(ncon) {
        if(mxGetN(pA) != ndec)
            mexErrMsgTxt("A has incompatible dimensions");
        if(mxGetNumberOfElements(pB) != ncon)
            mexErrMsgTxt("b has incompatible dimensions");
    }
    if(!mxIsEmpty(pLB) && (mxGetNumberOfElements(pLB) != ndec))
        mexErrMsgTxt("lb has incompatible dimensions");
    if(!mxIsEmpty(pUB) && (mxGetNumberOfElements(pUB) != ndec))
        mexErrMsgTxt("ub has incompatible dimensions");    
    if(nrhs > eXTYPE && !mxIsEmpty(pXTYPE) && (mxGetNumberOfElements(pXTYPE) != ndec))
        mexErrMsgTxt("xtype has incompatible dimensions");    
}

//Option Getting Methods
void getIntOption(const mxArray *opts, const char *option, int &var)
{
    if(mxGetField(opts,0,option))
        var = (int)*mxGetPr(mxGetField(opts,0,option));
}
void getDblOption(const mxArray *opts, const char *option, double &var)
{
    if(mxGetField(opts,0,option))
        var = *mxGetPr(mxGetField(opts,0,option));
}
void getStrOption(const mxArray *opts, const char *option, char *str)
{
    if(mxGetField(opts,0,option))
        mxGetString(mxGetField(opts,0,option),str,BUFSIZE);
}

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" SCIP-SDP: Solving Constraint Integer Programs Semidefinite Program Extension [v%s, Built %s]\n",SCIPSDP_VERSION,__DATE__);
    mexPrintf("  - Released under the Lesser GNU Public License\n");
    mexPrintf("  - Source available from: http://www.opt.tu-darmstadt.de/~smars/scip_sdp.html/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    mexPrintf("  - SCIP [%d.%d.%d] (ZIB Academic License)\n",SCIPmajorVersion(),SCIPminorVersion(),SCIPtechVersion());
    mexPrintf("  - SoPlex [v%d] (ZIB Academic License)\n",SOPLEX_VERSION);
    mexPrintf("  - DSDP [v%s] (Copyright 2004 University of Chicago)\n",DSDP_VERSION);

    mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}
