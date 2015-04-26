/* DSDPMEX - A MATLAB MEX Interface to DSDP
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.i2c2.aut.ac.nz
 */

/* Based in parts on dsdpmex.c supplied with DSDP 5.8 */

#include "mex.h"
#include "dsdp5.h"
#include "time.h"
#include "opti_util.h"

#ifdef LINK_MKL
    #include "mkl.h"
#endif

//Enable for Debug print out
// #define DEBUG

//DSDP Version
#define DSDP_VERSION "5.8"

//Macros & Defines
#define DSDP_ERR(x,msg) if(x) {sprintf(msgbuf,"%s, Error Code: %d",msg,x); mexErrMsgTxt(msgbuf);}
#define DSDP_MAX_TIME -27

//Argument Enums (in expected order of arguments)
enum {eF, eA, eB, eLB, eUB, eSDP, eY0, eOPTS};                   
//PRHS Defines    
#define pF      prhs[eF]
#define pA      prhs[eA]
#define pB      prhs[eB]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pSDP    prhs[eSDP]
#define pY0     prhs[eY0]
#define pOPTS   prhs[eOPTS]

//Function Prototypes
void printSolverInfo();
int addSDPCone(SDPCone sdpcone, const mxArray *cone, int block, double *X);
void checkInputs(const mxArray *prhs[], int nrhs);
void checkCone(const mxArray *cone, size_t ndec, size_t block);
void GetIntegerOption(const mxArray *opts, char *name, int *var);
void GetDoubleOption(const mxArray *opts, char *name, double *var);
static int DSDPMonitor(DSDP dsdp, void* dummy);

//Missing prototypes
extern int LPConeSetXVec(LPCone lpcone, double *xout, int n);

//Ctrl-C Detection
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" void utSetInterruptPending(bool);
#else
    extern bool utIsInterruptPending();
    extern void utSetInterruptPending(bool);
#endif
    
//Global Variables
char msgbuf[1024];
int printLevel;
double maxtime;
clock_t start, end;
    
//Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *f, *A = NULL, *b = NULL, *lb = NULL, *ub = NULL, *y0 = NULL;
    double *sdpDIM = NULL, *SDP_pr = NULL;
    
    //Return Args
    double *x, *pval, *dval, *exitflag, *iter, *pdflag;
    
    //Options (most get defaults written in)
    int maxiter = 1500;  
    int reuse=4,rpos=0,drho=1,ndim,sdpnmax=1;
    double penalty,rho,dbound,dlbound,zbar,r0,mu0,ylow,yhigh,gaptol,pnormtol,maxtrust,steptol,inftol,infptol;
    double lpb=1.0, datanorm[3], *dreuse, *fixed = NULL;
    
    //Internal Vars
    size_t nlincon = 0, ndec = 0, ncones = 0, nfix = 0;
    size_t lincon_nz = 0;    
    size_t i, j;
    size_t nLB = 0, nUB = 0;
    int *temp_ir = NULL, *temp_jc = NULL;
    double *temp_pr = NULL;
    const char *onames[2] = {"pval","dval"};
    const char *fnames[11] = {"iter","pdflag","r","mu","pstep","dstep","pnorm","ynorm","tracex","reuse","rho"};    
    double evaltime, *X = NULL;
    int iters = 0, status, indcell = 0;
    
    //DSDP Vars
    DSDP dsdp;
    SDPCone sdpcone = NULL;
    LPCone lpcone   = NULL;
    BCone bcone     = NULL;
    DSDPTerminationReason reason;
    DSDPSolutionType pdfeasible;  

    //Sparse Indicing
    mwIndex *A_ir, *A_jc;
    //Version Return
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
            plhs[0] = mxCreateString(DSDP_VERSION);
        return;
    }        
    
    //Check Inputs
    checkInputs(prhs,nrhs); 
    
    //Get pointers to Input variables
	f = mxGetPr(pF); ndec = mxGetNumberOfElements(pF);
    if(!mxIsEmpty(pA)) {
        A = mxGetPr(pA); 
        A_ir = mxGetIr(pA);
        A_jc = mxGetJc(pA);
        b = mxGetPr(pB);
        nlincon = mxGetM(pA);
        lincon_nz = A_jc[mxGetN(pA)];
    }
    if(nrhs > eLB && !mxIsEmpty(pLB))
        lb = mxGetPr(pLB); 
    if(nrhs > eUB && !mxIsEmpty(pUB))
        ub = mxGetPr(pUB);
    if(nrhs > eSDP && !mxIsEmpty(pSDP)) {
        if(mxIsCell(pSDP))
            ncones = mxGetNumberOfElements(pSDP);
        else
            ncones = 1;
    }
    if(nrhs > eY0 && !mxIsEmpty(pY0))
        y0 = mxGetPr(pY0);
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && mxGetField(pOPTS,0,"fixed") && !mxIsEmpty(mxGetField(pOPTS,0,"fixed"))) {
        fixed = mxGetPr(mxGetField(pOPTS,0,"fixed"));
        nfix = mxGetM(mxGetField(pOPTS,0,"fixed")); 
    }
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateStructMatrix(1,1,2,onames);
    mxSetField(plhs[1],0,onames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[1],0,onames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);   
    x = mxGetPr(plhs[0]); 
    pval = mxGetPr(mxGetField(plhs[1],0,onames[0]));
    dval = mxGetPr(mxGetField(plhs[1],0,onames[1]));
    exitflag = mxGetPr(plhs[2]);    
    //Info Output    
    plhs[3] = mxCreateStructMatrix(1,1,11,fnames);
    mxSetField(plhs[3],0,fnames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[2],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[3],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[4],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[5],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[6],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[7],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[8],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[9],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[10],mxCreateDoubleMatrix(1,1, mxREAL));
    iter = mxGetPr(mxGetField(plhs[3],0,fnames[0]));
    pdflag = mxGetPr(mxGetField(plhs[3],0,fnames[1]));   
    dreuse = mxGetPr(mxGetField(plhs[3],0,"reuse"));
    if(nlhs > 4)         
    	plhs[4] = mxCreateCellMatrix(ncones+(int)(nlincon>0)+(int)(nfix>0),1);        
    
    //Set Defaults
    maxtime = 1000;
    printLevel = 0;
    
    //Create DSDP Problem
    DSDP_ERR( DSDPCreate((int)ndec,&dsdp), "Error Creating DSDP Problem");
    //Set Monitor
    DSDP_ERR( DSDPSetMonitor(dsdp,DSDPMonitor,0), "Error Setting DSDP Monitor");
    
    //Set Dual Objective
    for (i=0;i<ndec;i++){
        DSDP_ERR( DSDPSetDualObjective(dsdp,(int)i+1,f[i]), "Error Adding Objective Coefficients"); }
    
    //Check finite bounds for allocation
    if(lb || ub)
        for(i=0;i<ndec;i++) {
            if(lb)
                if(!mxIsInf(lb[i]))
                    nLB++;
            if(ub)
                if(!mxIsInf(ub[i]))
                    nUB++;
        }
    
    //Set Bounds as BCone
    if(nLB || nUB) {
        DSDP_ERR( DSDPCreateBCone(dsdp, &bcone), "Error creating BCone");
        DSDP_ERR( BConeAllocateBounds(bcone, (int)(nLB+nUB)), "Error allocating bounds");        
        for(i=0;i<ndec;i++) {
            if(nLB > 0 && !mxIsInf(lb[i]))
                DSDP_ERR( BConeSetLowerBound(bcone, (int)i+1, lb[i]), "Error setting lower bound");
            if(nUB > 0 && !mxIsInf(ub[i]))
                DSDP_ERR( BConeSetUpperBound(bcone, (int)i+1, ub[i]), "Error setting upper bound");
        }
    }
    
    //Set Linear Inequality Constraints as LPCone
    if(nlincon) {
        int M = (int)mxGetM(pA);
        int N = (int)mxGetN(pA);
        DSDP_ERR( DSDPCreateLPCone(dsdp, &lpcone), "Error creating LPCone (inequalities)");
        //Create Memory to store A*x <= b in dsdp and integer format
        temp_jc = mxCalloc(N+2,sizeof(int));
        temp_ir = mxCalloc(lincon_nz+M,sizeof(int));
        temp_pr = mxCalloc(lincon_nz+M,sizeof(double));
        //Copy over linear A
        for(i=0;i<=(size_t)N;i++)
            temp_jc[i] = (int)A_jc[i];
        for(i=0;i<lincon_nz;i++) {
            temp_ir[i] = (int)A_ir[i];
            temp_pr[i] = A[i];
        }
        //Append linear rhs (b)
        temp_jc[N+1] = temp_jc[N] + M;
        for(i=lincon_nz,j=0;j<(size_t)M;j++) {
            if(b[j] != 0) {
                temp_ir[i] = (int)j;
                temp_pr[i++] = b[j];
            }
            else
                temp_jc[N+1]--;
        }
        #ifdef DEBUG
            mexPrintf("---- Inequality Constraints ----\n");
            for(i=0;i<=(size_t)(N+1);i++)
                mexPrintf("jc[%d] = %d\n",i,temp_jc[i]);
            for(i=0;i<lincon_nz+M;i++)
                mexPrintf("ir[%d] = %d, pr[%d] = %f\n",i,temp_ir[i],i,temp_pr[i]);
        #endif        
        //Set LP Cone Data
        DSDP_ERR( LPConeSetData2(lpcone, M, temp_jc, temp_ir, temp_pr), "Error setting LP Cone data (inequality)" );
        //Optionally set X data
        if(nlhs > 4) {
            mxSetCell(plhs[4],indcell,mxCreateDoubleMatrix(M,1,mxREAL));
            DSDP_ERR( LPConeSetXVec(lpcone,mxGetPr(mxGetCell(plhs[4],indcell++)),M), "Error setting LP Cone X data" );
        }
    }
    
    //Set Semidefinite Constraints as SDPCone
    if(ncones) {
        //Create the cone structure, specifying each constraint as a block
        DSDP_ERR( DSDPCreateSDPCone(dsdp,(int)ncones,&sdpcone), "Error creating SDPCone");
        //Add each constraint cone
        for(i=0;i<ncones;i++) {
            if(ncones == 1 && !mxIsCell(pSDP)) {
                if(nlhs > 4) {
                    mxSetCell(plhs[4],indcell,mxCreateDoubleMatrix(mxGetM(pSDP),1,mxREAL));
                    X = mxGetPr(mxGetCell(plhs[4],indcell++));
                }
                ndim = addSDPCone(sdpcone,pSDP,(int)i,X);
            }
            else {
                if(nlhs > 4) {
                    mxSetCell(plhs[4],indcell,mxCreateDoubleMatrix(mxGetM(mxGetCell(pSDP,i)),1,mxREAL));
                    X = mxGetPr(mxGetCell(plhs[4],indcell++));
                }
                ndim = addSDPCone(sdpcone,mxGetCell(pSDP,i),(int)i,X);
            }
            //Update max dim
            if(sdpnmax < ndim)
                sdpnmax = ndim;
        }
    }
    
    //Set y0
    if (y0)
        for (i=0;i<ndec;i++) {
            DSDP_ERR( DSDPSetY0(dsdp,(int)i+1,y0[i]), "Error setting Y0");            
        }
    
    //Determine whether to reuse schur complement matrix (dsdp authors' heuristic)
    if(ndec == 1)
        reuse = 1/sdpnmax;
    else
        reuse = ((int)ndec-2)/sdpnmax; 
    if (ndec<50 && reuse==0) reuse=1;
    if (reuse>=1) reuse++;
    reuse=reuse*reuse;
    if (ndec<2000 && ndec>10) reuse=10;
    if (ndec>12) reuse=12;    
    
    //Get DSDP Default Options
    DSDP_ERR( DSDPGetR(dsdp,&r0), "Error Getting R");
    DSDP_ERR( DSDPGetPenaltyParameter(dsdp,&penalty), "Error Getting Penalty Parameter");
    DSDP_ERR( DSDPGetPotentialParameter(dsdp,&rho), "Error Getting Potential Parameter");
    DSDP_ERR( DSDPGetDualBound(dsdp,&dbound), "Error Getting Dual Bound");
    DSDP_ERR( DSDPGetGapTolerance(dsdp,&gaptol), "Error Getting Gap Tolerance");
    DSDP_ERR( DSDPGetRTolerance(dsdp,&inftol), "Error Getting R Tolerance");
    DSDP_ERR( DSDPGetBarrierParameter(dsdp,&mu0), "Error Getting Barrier Parameter");
    DSDP_ERR( DSDPGetMaxTrustRadius(dsdp,&maxtrust), "Error Getting Max Trust Radius");
    DSDP_ERR( DSDPGetStepTolerance(dsdp,&steptol), "Error Getting Step Tolerance");
    DSDP_ERR( DSDPGetPTolerance(dsdp,&infptol), "Error Getting P Tolerance");
    DSDP_ERR( DSDPGetPNormTolerance(dsdp,&pnormtol), "Error Getting PNorm Tolerance");
    
    //Get Data Norms to establish y bounds
    DSDP_ERR( DSDPGetDataNorms(dsdp, datanorm), "Error Getting Data Norms");
    DSDP_ERR( DSDPGetYBounds(dsdp,&ylow,&yhigh), "Error Getting Y Bounds");
    if (datanorm[0]==0){DSDP_ERR( DSDPSetYBounds(dsdp,-1.0,1.0), "Error Setting Y Bounds");}
    
    //Get User Options (overwrites defaults above)
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        //OPTI Options
        GetIntegerOption(pOPTS,"maxiter",&maxiter);
        GetDoubleOption(pOPTS,"maxtime",&maxtime);
        GetIntegerOption(pOPTS,"display",&printLevel);
        //DSDP Options
        GetDoubleOption(pOPTS,"r0",&r0);
        GetDoubleOption(pOPTS,"penalty",&penalty);
        GetDoubleOption(pOPTS,"rho",&rho);
        GetDoubleOption(pOPTS,"dbound",&dbound);
        GetDoubleOption(pOPTS,"gaptol",&gaptol);
        GetDoubleOption(pOPTS,"rtol",&inftol);
        GetDoubleOption(pOPTS,"mu0",&mu0);
        GetDoubleOption(pOPTS,"maxtrust",&maxtrust);
        GetDoubleOption(pOPTS,"steptol",&steptol);
        GetDoubleOption(pOPTS,"ptol",&infptol);
        GetDoubleOption(pOPTS,"pnormtol",&pnormtol); 
        GetIntegerOption(pOPTS,"reuse",&reuse);
        GetIntegerOption(pOPTS,"rpos",&rpos);
        GetIntegerOption(pOPTS,"drho",&drho);
        //Check and set DSDP options without valid defaults
        if(mxGetField(pOPTS,0,"zbar") && !mxIsEmpty(mxGetField(pOPTS,0,"zbar"))) {
            GetDoubleOption(pOPTS,"zbar",&zbar);
            DSDP_ERR( DSDPSetZBar(dsdp,zbar), "Error Setting Z Bar");
        }
        if(mxGetField(pOPTS,0,"dlbound") && !mxIsEmpty(mxGetField(pOPTS,0,"dlbound"))) {
            GetDoubleOption(pOPTS,"dlbound",&dlbound);
            DSDP_ERR( DSDPSetDualLowerBound(dsdp,dlbound), "Error Setting Dual Lower Bound");
        }
        if(mxGetField(pOPTS,0,"ybound") && !mxIsEmpty(mxGetField(pOPTS,0,"ybound"))) {
            GetDoubleOption(pOPTS,"ybound",&yhigh); ylow = -yhigh;
            DSDP_ERR( DSDPSetYBounds(dsdp,ylow,yhigh), "Error Setting Y Bounds");
        }
    }

    //Set DSDP Options with Defaults
    DSDP_ERR( DSDPSetMaxIts(dsdp,maxiter), "Error Setting Max Iterations");    
    DSDP_ERR( DSDPSetR0(dsdp,r0), "Error Setting Option R0 ");        
    DSDP_ERR( DSDPSetPenaltyParameter(dsdp,penalty), "Error Setting Penalty Parameter");
    DSDP_ERR( DSDPSetPotentialParameter(dsdp,rho), "Error Setting Potential Parameter");
    DSDP_ERR( DSDPSetDualBound(dsdp,dbound), "Error Setting Dual Bound");
    DSDP_ERR( DSDPSetGapTolerance(dsdp,gaptol), "Error Setting Gap Tolerance");
    DSDP_ERR( DSDPSetRTolerance(dsdp,inftol), "Error Setting R Tolerance");
    DSDP_ERR( DSDPSetBarrierParameter(dsdp,mu0), "Error Setting Barrier Parameter");
	DSDP_ERR( DSDPSetMaxTrustRadius(dsdp,maxtrust), "Error Setting Max Trust Radius");
	DSDP_ERR( DSDPSetStepTolerance(dsdp,steptol), "Error Setting Step Tolerance")
    DSDP_ERR( DSDPSetPTolerance(dsdp,infptol), "Error Setting P Tolerance");
    DSDP_ERR( DSDPSetPNormTolerance(dsdp,pnormtol), "Error Setting PNorm Tolerance");   
    if(reuse < 0) reuse = 0; if(reuse > 15) reuse = 15;
    DSDP_ERR( DSDPReuseMatrix(dsdp,reuse), "Error Setting Reuse Matrix");    
    //Set Other DSDP Options
    DSDP_ERR( DSDPUsePenalty(dsdp,rpos), "Error Setting Use Penalty");
    DSDP_ERR( DSDPUseDynamicRho(dsdp,drho), "Error Setting Dynamic Rho");    
    if (lpb<0.1) lpb=0.1;
    if(lpcone) DSDP_ERR( LPConeScaleBarrier(lpcone,lpb), "Error Setting LPCone Scale Barrier");   
    
    //Set Fixed Variables
    if(fixed != NULL) {
        if(nlhs > 4) {            
            mxSetCell(plhs[4],indcell,mxCreateDoubleMatrix(nfix,1,mxREAL));
            X = mxGetPr(mxGetCell(plhs[4],indcell++));
        }
        else
            X = NULL;
        DSDP_ERR( DSDPSetFixedVariables(dsdp, fixed, &fixed[nfix], X, (int)nfix), "Error Setting Fixed Variables");
    }
    
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is DSDP v%s\n",DSDP_VERSION); 
        mexPrintf(" Authors: Steve Benson, Yinyu Ye and Xiong Zhang\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:        %4d\n",ndec);
        mexPrintf(" # Linear Inequalities:       %4d ",nlincon);
        if(nlincon)
            mexPrintf("[%d nz]\n",lincon_nz);
        else
            mexPrintf("\n");        
        mexPrintf(" # Semidefinite Cones:        %4d\n",ncones);

        mexPrintf("------------------------------------------------------------------\n");
    }
    
    //Start timer
    start = clock();
    //Call DSDP Setup to initialize problem
    DSDP_ERR( DSDPSetup(dsdp), "Error setting up DSDP Problem, likely out of memory");
    //Now Solve the Problem
    status = DSDPSolve(dsdp);
    //Stop Timer
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    //Determine Stop Reason
    if(status == 0) {
        DSDP_ERR( DSDPStopReason(dsdp,&reason), "Error retrieving post-solve stop reason"); }
    else if(status == DSDP_MAX_TIME || status == DSDP_USER_TERMINATION)
        reason = status;
    else {
        DSDP_ERR( status, "Error solving DSDP Problem!");}
    
    //Computer X and Get Solution Type
    if (reason!=DSDP_INFEASIBLE_START)
        DSDP_ERR( DSDPComputeX(dsdp), "Error computing post-solve x");
    DSDP_ERR( DSDPGetSolutionType(dsdp,&pdfeasible), "Error collecting post-solve solution type");
    
    //Copy Dual Solution
    DSDP_ERR( DSDPGetY(dsdp,x,(int)ndec), "Error returning Solution Vector");
    //Collect Output Statistics
    DSDPGetIts(dsdp,&iters);     
    DSDPGetDObjective(dsdp,dval);    
    DSDPGetPObjective(dsdp,pval);    
    DSDPGetR(dsdp,mxGetPr(mxGetField(plhs[3],0,"r")));
    DSDPGetBarrierParameter(dsdp,mxGetPr(mxGetField(plhs[3],0,"mu")));
    DSDPGetStepLengths(dsdp,mxGetPr(mxGetField(plhs[3],0,"pstep")),mxGetPr(mxGetField(plhs[3],0,"dstep")));
    DSDPGetPnorm(dsdp,mxGetPr(mxGetField(plhs[3],0,"pnorm")));
    DSDPGetYMaxNorm(dsdp,mxGetPr(mxGetField(plhs[3],0,"ynorm")));
    DSDPGetTraceX(dsdp,mxGetPr(mxGetField(plhs[3],0,"tracex")));
    DSDPGetPotentialParameter(dsdp,mxGetPr(mxGetField(plhs[3],0,"rho")));
    *dreuse = (double)reuse;
    
    //Assign to MATLAB
    *iter = (double)iters;    
    *exitflag = (double)reason;
    *pdflag = (double)pdfeasible;
    
    //Print Header
    if(printLevel){            
        //Detail termination reason
        switch(reason)
        {
            //Success
            case DSDP_CONVERGED:	
                mexPrintf("\n *** DSDP CONVERGED ***\n");  break;
            case DSDP_UPPERBOUND:
                mexPrintf("\n *** DSDP CONVERGED: Dual Objective exceeds its bound***\n"); break;
            //Error
            case DSDP_SMALL_STEPS:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Terminated due to Small Steps ***\n"); break;
            case DSDP_MAX_IT:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Maximum Iterations Reached ***\n"); break;
            case DSDP_MAX_TIME:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Maximum Time Reached ***\n"); break;
            case DSDP_INFEASIBLE_START:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Infeasible Starting Point ***\n"); break;
            case DSDP_USER_TERMINATION:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: User Exited ***\n"); break;            
            //Here is ok too?
            default:
                mexPrintf("\n *** DSDP FINISHED ***\n"); break;
        }
        //Detail solution status
        if(reason == DSDP_CONVERGED || reason == DSDP_UPPERBOUND) {
            switch(pdfeasible)
            {
                //Success
                case DSDP_PDFEASIBLE:
                    mexPrintf(" Solution Status: Both Primal and Dual are Feasible and Bounded\n"); break;             
                //Error
                case DSDP_UNBOUNDED:
                    mexPrintf(" Solution Status: Dual Unbounded, Primal Infeasible\n"); break;
                case DSDP_INFEASIBLE:
                    mexPrintf(" Solution Status: Primal Unbounded, Dual Infeasible\n"); break;
                case DSDP_PDUNKNOWN:
                default:
                    mexPrintf(" Solution Status: Unknown - Check Dual Bounds\n"); break;
            }
        }

        if(reason==DSDP_CONVERGED)
        	mexPrintf("\n Final Primal Objective:  %2.5g\n Final Dual Objective:    %2.5g\n In %5d iterations\n    %5.2f seconds\n",*pval,*dval,iters,evaltime);

        mexPrintf("------------------------------------------------------------------\n\n");
    }  

    //Free DSDP Problem
    DSDP_ERR( DSDPDestroy(dsdp), "Error Destroying DSDP Problem");
    //Free Temporary memory
    if(temp_jc)  {mxFree(temp_jc);  temp_jc  = NULL;}
    if(temp_ir)  {mxFree(temp_ir);  temp_ir  = NULL;}
    if(temp_pr)  {mxFree(temp_pr);  temp_pr  = NULL;}
}               

//Add SDP Constraint
int addSDPCone(SDPCone sdpcone, const mxArray *cone, int block, double *X)
{
    int i;
    double *SDP_pr  = mxGetPr(cone);
    mwIndex *SDP_ir = mxGetIr(cone);    
    mwIndex *SDP_jc = mxGetJc(cone);
    int *iSDP_ir    = NULL;
    int SDP_M       = (int)mxGetM(cone);
    int SDP_N       = (int)mxGetN(cone); //remember [C A0 A1 A2...] so non-square    
    int SDP_C_nnz   = 0; //nnz in C
    int SDP_A_nnz   = 0; //nnz in current A 
    int SDP_A_index = 0; //start index of pr and ir within current A
    int SDP_DIM     = (int)(sqrt(2*(double)SDP_M + 0.25) - 0.5); //calculate dimension based on rows (lower tri)
    
    #if DEBUG
        mexPrintf("SDP_DIM: %d, M: %d, N: %d\n",SDP_DIM,SDP_M,SDP_N);
    #endif
    //Copy cast ir
    iSDP_ir = (int*)mxCalloc(SDP_jc[SDP_N],sizeof(int));
    for(i=0;i<(int)SDP_jc[SDP_N];i++)
        iSDP_ir[i] = (int)SDP_ir[i];
        
    //Set Cone Size
    DSDP_ERR( SDPConeSetBlockSize(sdpcone,block,(int)SDP_DIM), "Error setting cone dimensions"); 
    DSDP_ERR( SDPConeUsePackedFormat(sdpcone, block), "Error setting cone packed format");
    DSDP_ERR( SDPConeSetSparsity(sdpcone,block,(int)SDP_jc[SDP_N]), "Error setting cone sparsity nnz");
    //Add C
    SDP_C_nnz = (int)(SDP_jc[1]-SDP_jc[0]); 
    #if DEBUG
        mexPrintf("C nnz: %d\n",SDP_C_nnz);
    #endif
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
    //Add X Vector if supplied
    if(X != NULL)
        DSDP_ERR( SDPConeSetXArray(sdpcone,block, SDP_DIM, X, SDP_M), "Error Setting SDPCone X Array" );
    return SDP_DIM;
    //Do not free memory here, MATLAB will take care of it once MEX completes (DSDP will use it during solving)  [See MATLAB/User's Guide/C/C++ and Fortran API Reference/mxFree] 
}

//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, m, i;
    const mxArray *fxd;
    double *val;
    
    //Correct number of inputs
    if(nrhs < 3)
        mexErrMsgTxt("You must supply at least 3 arguments to dsdp (f, A, b)"); 
    
    //Check we have an objective
    if(mxIsEmpty(pF))
        mexErrMsgTxt("You must supply an objective function via f!");
    
    //Check we have some constraints
    if(nrhs >= (eSDP+1) && mxIsEmpty(pA) && mxIsEmpty(pLB) && mxIsEmpty(pUB) && mxIsEmpty(pSDP))
        mexErrMsgTxt("You must supply constraints to this solver!");
    else if(nrhs == (eUB+1) && mxIsEmpty(pA) && mxIsEmpty(pLB) && mxIsEmpty(pUB))
        mexErrMsgTxt("You must supply constraints to this solver!");
    else if(nrhs == (eB+1) && (mxIsEmpty(pA) || mxIsEmpty(pB)))
        mexErrMsgTxt("You must supply constraints to this solver!");
    
    //Check options is a structure
    if(nrhs > eOPTS && !mxIsStruct(pOPTS))
        mexErrMsgTxt("The options argument must be a structure");
    
    //Get Sizes
    ndec = mxGetNumberOfElements(pF);
    
    //Check linear constraints
    if(!mxIsEmpty(pA)) {
        //Check types
        if(!mxIsSparse(pA) || mxIsComplex(pA) || !mxIsDouble(pA))
            mexErrMsgTxt("A must be real, sparse, double matrix");
        if(mxIsSparse(pB) || mxIsComplex(pB) || !mxIsDouble(pB))
            mexErrMsgTxt("b must be real, dense, double column vector");
        //Check sizes
        if(mxGetN(pA) != ndec)
            mexErrMsgTxt("A does not have the same number of columns as decision variables!");
        if(mxGetM(pA) != mxGetNumberOfElements(pB))
            mexErrMsgTxt("A and b do not have the same number of rows!");        
    }
    //Check bounds
    if(nrhs > eLB && !mxIsEmpty(pLB)) {
        if(mxIsSparse(pLB) || mxIsComplex(pLB) || !mxIsDouble(pLB))
            mexErrMsgTxt("lb must be real, dense, double column vector");
        if(mxGetNumberOfElements(pLB) != ndec)
            mexErrMsgTxt("lb is not the same length as f!");
    }
    if(nrhs > eUB && !mxIsEmpty(pUB)) {
        if(mxIsSparse(pUB) || mxIsComplex(pUB) || !mxIsDouble(pUB))
            mexErrMsgTxt("ub must be real, dense, double column vector");
        if(mxGetNumberOfElements(pUB) != ndec)
            mexErrMsgTxt("ub is not the same length as f!");
    }
    //Check sdcones
    if(nrhs > eSDP && !mxIsEmpty(pSDP)) {
        if(mxIsCell(pSDP)) {
            for(i=0;i<mxGetNumberOfElements(pSDP);i++)
                checkCone(mxGetCell(pSDP,i),ndec,i);
        }
        else
            checkCone(pSDP,ndec,1);        
    }
    //Check y0
    if(nrhs > eY0 && !mxIsEmpty(pY0)) {
        if(mxIsSparse(pY0) || mxIsComplex(pY0) || !mxIsDouble(pY0))
            mexErrMsgTxt("y0 must be real, dense, double column vector");
        if(mxGetNumberOfElements(pY0) != ndec)
            mexErrMsgTxt("y0 is not the same length as f!");
    }
    //Check fixed vars option
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && mxGetField(pOPTS,0,"fixed") && !mxIsEmpty(mxGetField(pOPTS,0,"fixed"))) {
        fxd = mxGetField(pOPTS,0,"fixed");
        //Check Type
        if(mxIsSparse(fxd) || mxIsComplex(fxd) || !mxIsDouble(fxd))
            mexErrMsgTxt("The fixed option must be a real, dense, double matrix");
        val = mxGetPr(fxd);
        m = mxGetM(fxd);
        if(mxGetN(fxd) != 2)
            mexErrMsgTxt("The fixed option must contain 2 columns");
        if(m > ndec)
            mexErrMsgTxt("The fixed option must not contain more rows than decision variables");
        for(i=0;i<m;i++) {
            if(val[i] <= 0 || val[i] > ndec)
                mexErrMsgTxt("A variable index in the fixed option is <= 0 or > ndec");
            if(mxIsNaN(val[i+m]) || mxIsInf(val[i+m]))
                mexErrMsgTxt("A variable value in the fixed option is NaN or Inf");            
        }
    }
}

//Check cone for correct dims and data type
void checkCone(const mxArray *cone, size_t ndec, size_t block)
{
    char msgbuf[1024];
    size_t M;
    double Msq;
    //Check type
    if(!mxIsSparse(cone) || mxIsComplex(cone) || !mxIsDouble(cone)) {
        sprintf(msgbuf,"Cone %d must be a real, sparse, double matrix",block);
        mexErrMsgTxt(msgbuf);
    }
    //Check we can make a square matrix from M
    M = mxGetM(cone);
    Msq = floor(sqrt(2.0*(double)M));
    if((Msq*(Msq+1)/2.0) != (double)M) {
        sprintf(msgbuf,"Cone %d does not contain a number of rows that can be converted to a square matrix! DSDPMEX expects just TRIU elements.",block);
        mexErrMsgTxt(msgbuf);
    }       
    //Check correct number of columns
    if(mxGetN(cone) != (ndec+1)) {
        sprintf(msgbuf,"Cone %d does not contain the correct number of columns! Expected %d",block,ndec+1);
        mexErrMsgTxt(msgbuf);
    }
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

//Solver Iteration Callback Monitor
static int DSDPMonitor(DSDP dsdp, void* dummy)
{
    int    iter;
    double pobj,dobj,pstp=0,dstp,mu,res,pnorm,pinfeas;
    double evaltime;
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
  
    //Check max time
    if(evaltime > maxtime)
        return DSDP_MAX_TIME;

    //Check ctrl-c
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting DSDP...\n\n");
        return DSDP_USER_TERMINATION; //terminate
    }

    if(printLevel>1) {
        DSDPGetIts(dsdp,&iter);
        DSDPGetDDObjective(dsdp,&dobj);
    	DSDPGetPPObjective(dsdp,&pobj);
    	DSDPGetR(dsdp,&res);
    	DSDPGetPInfeasibility(dsdp,&pinfeas);
    	DSDPGetStepLengths(dsdp,&pstp,&dstp);
    	DSDPGetBarrierParameter(dsdp,&mu);
    	DSDPGetPnorm(dsdp,&pnorm);
        
        //Display heading if % 20 iters
        if(iter==0 || iter % 20 == 0)
            mexPrintf("Iter   Time[s]   PP Objective      DD Objective     PInfeas   DInfeas     Nu     StepLength   Pnrm\n");
            
        //Display parameters
        mexPrintf("%-3d   %6.2f  %16.8e  %16.8e  %9.1e %9.1e %9.1e",iter,evaltime,pobj,dobj,pinfeas,res,mu);
        mexPrintf("  %4.2f  %4.2f",pstp,dstp);
        if (pnorm>1.0e3)
            mexPrintf("  %1.0e \n",pnorm);
        else
            mexPrintf("  %5.2f \n",pnorm);
        mexEvalString("drawnow;"); //flush draw buffer
    }
    //Return ok
    return(0);
}

//Print Solver Information
void printSolverInfo()
{    
    char vbuf[6]; getVSVer(vbuf);  
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" DSDP: Software for Semidefinite Programming [v%s, Built %s, VS%s]\n",DSDP_VERSION,__DATE__,vbuf);
    mexPrintf("  - Copyright 2004 University of Chicago: http://www.mcs.anl.gov/hs/software/DSDP/Copyright.txt\n");
    mexPrintf("  - Source available from: http://www.mcs.anl.gov/hs/software/DSDP/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    #endif

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.i2c2.aut.ac.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}