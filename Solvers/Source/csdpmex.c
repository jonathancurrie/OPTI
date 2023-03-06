/* CSDPMEX - A MATLAB MEX Interface to CSDP
 * Released Under the BSD 3-Clause License:
 * https://www.controlengineering.co.nz/Wikis/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2018
 * www.inverseproblem.co.nz
 */

/* NOTE: This interface expects the primal form SDP */

#include "mex.h"
#include "declarations.h"
#include "time.h"
#include "opti_build_utils.h"

#ifdef LINK_MKL
    #include "mkl.h"
#endif

//Enable for Debug print out
// #define DEBUG

//CSDP Defines
#define CSDP_VERSION "6.2 beta"
#define CSDP_MAX_TIME -27
#define CSDP_USER_TERMINATION -50
//Memory alignment macro
#define malloc(x) _aligned_malloc(x,16)

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
int addBounds(struct blockmatrix C, double *lb, double *ub, size_t ndec, int *nLB, int *nUB);
int insertBound(struct sparseblock *blockptr, double *bnd, int nBND, int con_num, int block, int ind, double coeff);
void insertLPVector(struct sparseblock *blockptr, const mxArray *Alin, int con_num, int block);
int addCMatrix(struct blockmatrix C, const mxArray *cone, int block);
void addAMatrix(struct sparseblock *blockptr, const mxArray *cone, int con_num, int block);
void checkInputs(const mxArray *prhs[], int nrhs);
void checkCone(const mxArray *cone, size_t ndec, size_t block);
void GetIntegerOption(const mxArray *opts, char *name, int *var);
void GetDoubleOption(const mxArray *opts, char *name, double *var);

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
int citer,printLevel;
double maxtime;
clock_t start, end;
    
//Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    //Input Args
    double *f, *blin = NULL, *lb = NULL, *ub = NULL, *y0 = NULL;
    
    //Return Args
    double *x, *pval, *dval, *exitflag, *iter, *pinf, *dinf, *realgap, *xzgap;
    
    //Options (most get defaults written in)
    int maxiter = 1500;  
    
    //Internal Vars
    size_t ndec = 0, nlincon = 0, lincon_nz = 0, total_dim = 0, ncones = 0; 
    size_t i, j;
    const char *onames[2] = {"pval","dval"};
    const char *fnames[5] = {"iter","pinf","dinf","realgap","xzgap"};    
    double evaltime;
    int status = -1, nb = 0, linoffset = 0, indlb = 1, indub = 1, nLB = 0, nUB = 0;
    mwIndex *jc;
    
    //CSDP Problem Data
    struct blockmatrix C;
    double *b, *y, *xx, objconstant = 0.0;
    struct constraintmatrix *constraints;
    struct blockmatrix X, Z;
    struct sparseblock *blockptr;
    struct paramstruc params;

    //Version Return
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString(CSDP_VERSION);
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }        
    
    //Check Inputs
    checkInputs(prhs,nrhs); 
    
    //Get pointers to Input variables
	f = mxGetPr(pF); ndec = mxGetNumberOfElements(pF);
    if(!mxIsEmpty(pA)) {
        blin = mxGetPr(pB);
        nlincon = mxGetM(pA);
        jc = mxGetJc(pA);
        lincon_nz = jc[ndec];        
    }
    if(nrhs > eLB && !mxIsEmpty(pLB)) {
        lb = mxGetPr(pLB); 
        //Ensure we have at least one finite bound
        for(i=0,j=0;i<ndec;i++)
            if(mxIsInf(lb[i]))
                j++;
        if(j==ndec)
            lb = NULL;
    }
    if(nrhs > eUB && !mxIsEmpty(pUB)) {
        ub = mxGetPr(pUB);
        //Ensure we have at least one finite bound
        for(i=0,j=0;i<ndec;i++)
            if(mxIsInf(ub[i]))
                j++;
        if(j==ndec)
            ub = NULL;
    }
    if(nrhs > eSDP && !mxIsEmpty(pSDP)) {
        if(mxIsCell(pSDP))
            ncones = mxGetNumberOfElements(pSDP);
        else
            ncones = 1;
    }
    if(nrhs > eY0 && !mxIsEmpty(pY0))
        y0 = mxGetPr(pY0);
    
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
    plhs[3] = mxCreateStructMatrix(1,1,5,fnames);
    mxSetField(plhs[3],0,fnames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[2],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[3],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[4],mxCreateDoubleMatrix(1,1, mxREAL));
    iter = mxGetPr(mxGetField(plhs[3],0,fnames[0]));
    pinf = mxGetPr(mxGetField(plhs[3],0,fnames[1]));
    dinf = mxGetPr(mxGetField(plhs[3],0,fnames[2]));
    realgap = mxGetPr(mxGetField(plhs[3],0,fnames[3]));
    xzgap = mxGetPr(mxGetField(plhs[3],0,fnames[4]));
    
    //Set Defaults
    citer = 0;
    maxtime = 1000;
    printLevel = 0;
    
    //Allocate Initial Storage for the Problem
    b = (double*)malloc((ndec+1)*sizeof(double));       //objective vector
    //C matrices [LB UB LIN SD]
    nb = (int)ncones+(int)(nlincon>0)+(int)(lb!=NULL)+(int)(ub!=NULL);
    #ifdef DEBUG
        mexPrintf("Number of blocks (including bounds, linear and sdcones): %d\n",nb);
    #endif            
    C.nblocks   = nb;
    C.blocks    = (struct blockrec*)malloc((nb+1)*sizeof(struct blockrec)); //+1 due to fortran index
    if(C.blocks == NULL)
        mexErrMsgTxt("Error allocating memory for C matrices");
    //Constraints (i.e. 1 per decision variable)
    constraints = (struct constraintmatrix*)malloc((ndec+1)*sizeof(struct constraintmatrix)); //+1 due to fortran index
    if(constraints == NULL) {
        free(C.blocks);
        mexErrMsgTxt("Error allocating memory for A matrices");
    }
    for(i=1;i<=ndec;i++)
        constraints[i].blocks=NULL; //initially set as NULL
    
    //Copy in and negate objective vector
    for(i=0;i<ndec;i++)
        b[i+1] = -f[i];
    
    //Create Bounds if Present
    if(lb || ub)
        linoffset += addBounds(C,lb,ub,ndec,&nLB,&nUB);
    
    //Create Linear Cone (diagonal) if present
    if(nlincon) {
        //Initialize C
        C.blocks[linoffset+1].blockcategory = DIAG;
        C.blocks[linoffset+1].blocksize = (int)nlincon;
        C.blocks[linoffset+1].data.vec = (double*)malloc((nlincon+1)*sizeof(double));
        if(C.blocks[linoffset+1].data.vec == NULL)
            mexErrMsgTxt("Error allocating memory for LP C diagonal");
        #ifdef DEBUG
            mexPrintf("LP C[%d] Vector size: %d\n",linoffset+1,nlincon);
        #endif        
        //Copy Elements
        for(i=0;i<nlincon;i++) {
            C.blocks[linoffset+1].data.vec[i+1] = -blin[i];
            #ifdef DEBUG
                mexPrintf("  C vec[%d] = %f\n",i+1,-blin[i]);
            #endif
        }
        linoffset++;
    }
    
    #ifdef DEBUG
        mexPrintf("\nBlock offset after bounds + linear con: %d\n\n",linoffset);
    #endif
    
    //Setup Semidefinite C matrices (note all full matrices, dense, in order from 1)
    for(i=1;i<=ncones;i++) {
        //Single Cone
        if(ncones == 1 && !mxIsCell(pSDP))
            total_dim += addCMatrix(C,pSDP,(int)i+linoffset);
        //Multiple Cones
        else 
            total_dim += addCMatrix(C,mxGetCell(pSDP,i-1),(int)i+linoffset);     
    }    
    //Add Linear Dims
    total_dim += nLB+nUB+nlincon;
    
    #ifdef DEBUG
        mexPrintf("\nTotal dimension of all cones: %d\n\n",total_dim);
    #endif
    
    //Setup Each Constraint (for each decision var) (in order from 1)
    indlb = 1; indub = 1;
    for(i=1;i<=ndec;i++) { 
        //For each Semidefinte A matrix (sparse triu, in reverse order, i.e. [SD, LP, UB, LB])
        for(j=ncones;j>0;j--) {
            //Create an A matrix            
            blockptr=(struct sparseblock*)malloc(sizeof(struct sparseblock));
            if(blockptr==NULL) {
                sprintf(msgbuf,"Error allocating memory for Semidefinite A[%d,%d]",i,j+linoffset);
                mexErrMsgTxt(msgbuf);
            }            
            //Single Cone
            if(ncones == 1 && !mxIsCell(pSDP)) 
                addAMatrix(blockptr,pSDP,(int)i,(int)j+linoffset); 
            //Multiple Cones
            else 
                addAMatrix(blockptr,mxGetCell(pSDP,j-1),(int)i,(int)j+linoffset);
            
            //Insert A matrix into constraint list
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }  
        
        //Linear Inequality Constraints
        if(nlincon) {
            //Create an A matrix
             blockptr=(struct sparseblock*)malloc(sizeof(struct sparseblock));
            if(blockptr==NULL) {
                sprintf(msgbuf,"Error allocating memory for LP A[%d]",i);
                mexErrMsgTxt(msgbuf);
            }
            //Insert LP A entries
            j = 1 + (int)(nUB > 0) + (int)(nLB > 0);
            insertLPVector(blockptr, pA, (int)i, (int)j);             
            //Insert A matrix into constraint list
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }
        
        //Upper Bounds
        if(nUB) {
            //Create an A matrix
            blockptr=(struct sparseblock*)malloc(sizeof(struct sparseblock));
            if(blockptr==NULL) {
                sprintf(msgbuf,"Error allocating memory for UB A[%d]",i);
                mexErrMsgTxt(msgbuf);
            }            
            //Insert Bound A matrix entries
            if(nLB > 0)
                indub += insertBound(blockptr,ub,nUB,(int)i,2,indub,1.0); //block 2 ([LB,UB,..] 1.0 for ub)
            else
                indub += insertBound(blockptr,ub,nUB,(int)i,1,indub,1.0); //block 1 (first block, 1.0 for ub)
            //Insert A matrix into constraint list
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }
        
        //Lower Bounds
        if(nLB) {
            //Create an A matrix
            blockptr=(struct sparseblock*)malloc(sizeof(struct sparseblock));
            if(blockptr==NULL) {
                sprintf(msgbuf,"Error allocating memory for LB A[%d]",i);
                mexErrMsgTxt(msgbuf);
            }            
            //Insert Bound A matrix entries
            indlb += insertBound(blockptr,lb,nLB,(int)i,1,indlb,-1.0); //block 1 (always first block, -1.0 for lb)
            //Insert A matrix into constraint list
            blockptr->next=constraints[i].blocks;
            constraints[i].blocks=blockptr;
        }
    }

//     //Set y0
//     if (y0)
//         for (i=0;i<ndec;i++) {
//             DSDP_ERR( DSDPSetY0(dsdp,(int)i+1,y0[i]), "Error setting Y0");            
//         }
//     
    
    //Get CSDP Default Options
    initparams(&params,&printLevel); 
    //Set OPTI default printLevel (none)
    printLevel = 0;
    
    //Get User Options (overwrites defaults above)
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS)) {
        //OPTI Options
        GetIntegerOption(pOPTS,"maxiter",&params.maxiter);
        GetDoubleOption(pOPTS,"maxtime",&maxtime);
        GetIntegerOption(pOPTS,"display",&printLevel);
        GetDoubleOption(pOPTS,"objconstant",&objconstant);
        //CSDP Options
        GetDoubleOption(pOPTS,"axtol",&params.axtol);
        GetDoubleOption(pOPTS,"atytol",&params.atytol);
        GetDoubleOption(pOPTS,"objtol",&params.objtol);
        GetDoubleOption(pOPTS,"pinftol",&params.pinftol);
        GetDoubleOption(pOPTS,"dinftol",&params.dinftol);
        GetDoubleOption(pOPTS,"minstepfrac",&params.minstepfrac);
        GetDoubleOption(pOPTS,"maxstepfrac",&params.maxstepfrac);
        GetDoubleOption(pOPTS,"minstepp",&params.minstepp);
        GetDoubleOption(pOPTS,"minstepd",&params.minstepd);
        GetIntegerOption(pOPTS,"usexzgap",&params.usexzgap);
        GetIntegerOption(pOPTS,"tweakgap",&params.tweakgap); 
        GetIntegerOption(pOPTS,"affine",&params.affine);
        GetDoubleOption(pOPTS,"perturbobj",&params.perturbobj);
        CheckOptiVersion(pOPTS);
        //Optionally write problem to a SDPA sparse file
        if(mxGetField(pOPTS,0,"writeprob") && !mxIsEmpty(mxGetField(pOPTS,0,"writeprob")) && mxIsChar(mxGetField(pOPTS,0,"writeprob"))) {
            mxGetString(mxGetField(pOPTS,0,"writeprob"),msgbuf,1024);
            write_prob(msgbuf,(int)total_dim,(int)ndec,C,b,constraints);
        }
    }

    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        mexPrintf(" This is CSDP v%s\n",CSDP_VERSION); 
        mexPrintf(" Author: Brian Borchers\n\n");
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
    //Find Initial Solution
    initsoln((int)total_dim,(int)ndec,C,b,constraints,&X,&y,&Z);
    //Solve the problem
    status=easy_sdp((int)total_dim,(int)ndec,C,b,constraints,objconstant,params,&X,&y,&Z,pval,dval,pinf,dinf,realgap,xzgap);
    //Stop Timer
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Copy and negate solution
    for(i=0;i<ndec;i++)
        x[i] = -y[i+1];
    //Assign other MATLAB outputs
    *iter = (double)citer-1;    
    *exitflag = (double)status;
     
    //Print Header
    if(printLevel){            
        //Detail termination reason
        switch(status)
        {
            //Success
            case 0:	
                mexPrintf("\n *** CSDP CONVERGED ***\n");  break;
            //Infeasible
            case 1:
                mexPrintf("\n *** TERMINATION: Primal Infeasible ***\n");  break;
            case 2:
                mexPrintf("\n *** TERMINATION: Dual Infeasible ***\n");  break;
            //Partial Success
            case 3:
                mexPrintf("\n *** TERMINATION: PARTIAL SUCESS ***\n *** A Solution is found but full accuracy was not achieved ***\n"); break;
            //Error
            case 4:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Maximum Iterations Reached ***\n"); break;
            case CSDP_MAX_TIME:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Maximum Time Reached ***\n"); break;
            case 5:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Stuck at edge of primal feasibility ***\n"); break;
            case 6:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Stuck at edge of dual infeasibility ***\n"); break;
            case 7:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Lack of progress ***\n"); break;
            case 8:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: X, Z, or O was singular ***\n"); break;
            case 9:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Detected NaN or Inf values ***\n"); break;
            case 10:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Easy-SDP General Failure ***\n"); break;
            case 11:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Failed C Check - Check Symmetry! ***\n"); break;
            case 12:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Failed Constraints Check ***\n"); break;
            case CSDP_USER_TERMINATION:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: User Exited ***\n"); break;            
            //Here is ok too?
            default:
                mexPrintf("\n *** CSDP FINISHED ***\n"); break;
        }

        if(status==0 || status==3)
        	mexPrintf("\n Final Primal Objective:  %2.5g\n Final Dual Objective:    %2.5g\n In %5d iterations\n    %5.2f seconds\n",*pval,*dval,citer-1,evaltime);

        mexPrintf("------------------------------------------------------------------\n\n");
    }  
    
    //Optionally write solution to a SDPA sparse file
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && mxGetField(pOPTS,0,"writesol") && !mxIsEmpty(mxGetField(pOPTS,0,"writesol")) && mxIsChar(mxGetField(pOPTS,0,"writesol"))) {
        mxGetString(mxGetField(pOPTS,0,"writesol"),msgbuf,1024);
        write_sol(msgbuf,(int)total_dim,(int)ndec,X,y,Z);
    }
    
    //Optionally retrieve X
    if(nlhs > 4) {
        plhs[4] = mxCreateCellMatrix(X.nblocks,1);
        for(i=0;i<X.nblocks;i++) {            
            //Set Block Values
            if(X.blocks[i+1].blockcategory == DIAG) {
                mxSetCell(plhs[4],i,mxCreateDoubleMatrix(X.blocks[i+1].blocksize,1,mxREAL)); //create vector
                xx = mxGetPr(mxGetCell(plhs[4],i));
                for(j=0;j<X.blocks[i+1].blocksize;j++)
                    xx[j] = X.blocks[i+1].data.vec[j+1]; 
            }
            else {
                mxSetCell(plhs[4],i,mxCreateDoubleMatrix(X.blocks[i+1].blocksize,X.blocks[i+1].blocksize,mxREAL)); //create matrix
                xx = mxGetPr(mxGetCell(plhs[4],i));
                for(j=0;j<(X.blocks[i+1].blocksize*X.blocks[i+1].blocksize);j++)
                    xx[j] = X.blocks[i+1].data.mat[j];           
            }
        }        
    }

    //Free CSDP Problem (including all allocated memory)
    free_prob((int)total_dim,(int)ndec,C,b,constraints,X,y,Z);
}    

//Create Finite Lower and Upper Bounds
int addBounds(struct blockmatrix C, double *lb, double *ub, size_t ndec, int *nLB, int *nUB)
{
    size_t i, ind;
    int  nb = 0;
    *nLB = 0; *nUB = 0;
        
    //Check finite bounds for allocation
    for(i=0;i<ndec;i++) {
        if(lb)
            if(!mxIsInf(lb[i]))
                (*nLB)++;
        if(ub)
            if(!mxIsInf(ub[i]))
                (*nUB)++;
    }
        
    //Create C vector (just diagonal) for lb
    if(*nLB) {
        //Initialize C
        C.blocks[nb+1].blockcategory = DIAG;
        C.blocks[nb+1].blocksize = *nLB;
        C.blocks[nb+1].data.vec = (double*)malloc((*nLB+1)*sizeof(double));
        if(C.blocks[nb+1].data.vec == NULL)
            mexErrMsgTxt("Error allocating memory for LB C diagonal");
        #ifdef DEBUG
            mexPrintf("LB C[%d] Vector size: %d\n",nb+1,*nLB);
        #endif
        //Copy Elements
        ind=1;
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(lb[i])) {
                C.blocks[nb+1].data.vec[ind] = lb[i];
                #ifdef DEBUG
                    mexPrintf("  C vec[%d] = %f\n",ind,lb[i]);
                #endif
                ind++;
            }
        }
        nb++;
    }
    //Create C vector (just diagonal) for ub
    if(*nUB) {
        //Initialize C
        C.blocks[nb+1].blockcategory = DIAG;
        C.blocks[nb+1].blocksize = *nUB;
        C.blocks[nb+1].data.vec = (double*)malloc((*nUB+1)*sizeof(double));
        if(C.blocks[nb+1].data.vec == NULL)
            mexErrMsgTxt("Error allocating memory for UB C diagonal");
        #ifdef DEBUG
            mexPrintf("UB C[%d] Vector size: %d\n",nb+1,*nUB);
        #endif
        //Copy Elements
        ind=1;
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(ub[i])) {
                C.blocks[nb+1].data.vec[ind] = -ub[i];
                #ifdef DEBUG
                    mexPrintf("  C vec[%d] = %f\n",ind,-ub[i]);
                #endif
                ind++;
            }
        }
        nb++;
    }
    return nb;
}

//Insert Bound Vector
int insertBound(struct sparseblock *blockptr, double *bnd, int nBND, int con_num, int block, int ind, double coeff)
{  
    //Assign Sizes
    blockptr->blocknum = block;
    blockptr->blocksize = nBND;
    blockptr->constraintnum = con_num;
    blockptr->next = NULL;
    blockptr->nextbyblock = NULL;
    //Check we have a finite constraint on this variable
    if(!mxIsInf(bnd[con_num-1])) {
        blockptr->numentries = 1;
        //Allocate LB A Memory (+1s due to Fortran index)
        blockptr->entries   = (double*)malloc((1+1)*sizeof(double));
        if(blockptr->entries == NULL) {
            if(coeff==-1.0)
                sprintf(msgbuf,"Error allocating entries memory for LB[%d, block %d]",con_num, block);
            else
                sprintf(msgbuf,"Error allocating entries memory for UB[%d, block %d]",con_num, block);
            mexErrMsgTxt(msgbuf);
        }
        blockptr->iindices  = (unsigned*)malloc((1+1)*sizeof(unsigned));
        if(blockptr->iindices == NULL) {
            if(coeff==-1.0)
                sprintf(msgbuf,"Error allocating iindices memory for LB[%d, block %d]",con_num, block);
            else
                sprintf(msgbuf,"Error allocating iindices memory for UB[%d, block %d]",con_num, block);
            mexErrMsgTxt(msgbuf);
        }
        blockptr->jindices  = (unsigned*)malloc((1+1)*sizeof(unsigned));
        if(blockptr->jindices == NULL) {
            if(coeff==-1.0)
                sprintf(msgbuf,"Error allocating jindices memory for LB[%d, block %d]",con_num, block);
            else
                sprintf(msgbuf,"Error allocating jindices memory for UB[%d, block %d]",con_num, block);
            mexErrMsgTxt(msgbuf);
        }
        //Fill in entry
        blockptr->entries[1]  = coeff;
        blockptr->iindices[1] = ind;    //don't use con_num incase leading elements are inf (lb = [-Inf,0] .. etc)
        blockptr->jindices[1] = ind;
        #ifdef DEBUG
            if(coeff==-1.0)
                mexPrintf("LB[%d, block %d] i: %d, j: %d, val: %f\n",con_num,block,ind,ind,coeff);
            else
                mexPrintf("UB[%d, block %d] i: %d, j: %d, val: %f\n",con_num,block,ind,ind,coeff);
        #endif
    }
    //Else skip entry
    else {
        blockptr->numentries = 0;
        blockptr->entries  = NULL;
        blockptr->iindices = NULL;
        blockptr->jindices = NULL;
    }
    return blockptr->numentries;
}

//Add LP A Vector
void insertLPVector(struct sparseblock *blockptr, const mxArray *Alin, int con_num, int block)
{
    size_t i, M;
    int nz = 0, ind = 0;    
    mwIndex *jc, *ir;
    double *pr;
    
    //Get Sparse Alin Matrix from Matlab
    ir  = mxGetIr(Alin);
    jc  = mxGetJc(Alin);
    pr  = mxGetPr(Alin);
    M   = mxGetM(Alin); 

    //Assign Sizes to CSDP A Matrix
    blockptr->blocknum = block;
    blockptr->blocksize = (int)M;
    blockptr->constraintnum = con_num;
    blockptr->next = NULL;
    blockptr->nextbyblock = NULL;
    
    //Determine number of nz elements in this column A(con_num - 1)
    nz = (int)(jc[con_num]-jc[con_num-1]);
    blockptr->numentries = nz;
    
    //Allocate A Memory (+1s due to Fortran index)
    blockptr->entries   = (double*)malloc((nz+1)*sizeof(double));
    if(blockptr->entries == NULL) {
        sprintf(msgbuf,"Error allocating entries memory for LP A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    blockptr->iindices  = (unsigned*)malloc((nz+1)*sizeof(unsigned));
    if(blockptr->iindices == NULL) {
        sprintf(msgbuf,"Error allocating iindices memory for LP A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    blockptr->jindices  = (unsigned*)malloc((nz+1)*sizeof(unsigned));
    if(blockptr->jindices == NULL) {
        sprintf(msgbuf,"Error allocating jindices memory for LP A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    
    #ifdef DEBUG
        mexPrintf("LP A[%d, block %d] dim: %d, nz: %d: \n",con_num,block,M,nz);
    #endif 
        
    //Fill in A values
    ind = 1; //Fortran start
    for(i=jc[con_num-1];i<jc[con_num];i++) {
        blockptr->iindices[ind] = (int)ir[i]+1;
        blockptr->jindices[ind] = (int)ir[i]+1;
        blockptr->entries[ind] = pr[i];

        #ifdef DEBUG
            mexPrintf("  iind[%d] = %d, jind[%d] = %d, entr[%d] = %f\n",ind,ir[i]+1,ind,ir[i]+1,ind,pr[i]);
        #endif
        ind++;
    } 
}

//Add CSDP C Matrix [dense, full matrix]
int addCMatrix(struct blockmatrix C, const mxArray *cone, int block)
{
    size_t i, M;
    int dim;
    mwIndex *jc, *ir;
    double *pr;
    
    //Get Sparse C Matrix from Matlab
    ir  = mxGetIr(cone);
    jc  = mxGetJc(cone);
    pr  = mxGetPr(cone);
    M   = mxGetM(cone); 
    dim = (int)sqrt(M);
    //Assign to CSDP C
    C.blocks[block].blockcategory = MATRIX;
    C.blocks[block].blocksize = dim;         
    C.blocks[block].data.mat = (double*)malloc(M*sizeof(double));
    if(C.blocks[block].data.mat == NULL) {
        sprintf(msgbuf,"Error allocating memory for C[%d]",block);
        mexErrMsgTxt(msgbuf);
    }
    #ifdef DEBUG
        mexPrintf("SD C[%d]: dim %d, nnz %d\n",block,dim,jc[1]-jc[0]);
    #endif
    //Initialise all elements to zero
    for(i=0;i<M;i++)
        C.blocks[block].data.mat[i] = 0.0;    
    //Copy in nz Elements (just first column contains C)
    for(i=jc[0]; i < jc[1]; i++) {
        C.blocks[block].data.mat[ir[i]] = pr[i];
        #ifdef DEBUG
            mexPrintf("C.mat[%d] = %f\n",ir[i],pr[i]);
        #endif
    }
    return dim;
}

//Add CSDP A Matrix [sparse, triu matrix]
void addAMatrix(struct sparseblock *blockptr, const mxArray *cone, int con_num, int block)
{
    size_t i, M;
    int dim, triunz = 0, cind, rind, ind = 0;    
    mwIndex *jc, *ir;
    double *pr;
    
    //Get Sparse [C A0 A1 A2] Matrix from Matlab
    ir  = mxGetIr(cone);
    jc  = mxGetJc(cone);
    pr  = mxGetPr(cone);
    M   = mxGetM(cone); 
    dim = (int)sqrt(M);

    //Assign Sizes to CSDP A Matrix
    blockptr->blocknum = block;
    blockptr->blocksize = dim;
    blockptr->constraintnum = con_num;
    blockptr->next = NULL;
    blockptr->nextbyblock = NULL;
    
    //Determine number of nz elements in upper triangle
    for(i=jc[con_num];i<jc[con_num+1];i++) { //read down the nzs in the A column
        //Determine row index
        rind = ir[i] % dim;
        //If row <= col, in upper tri
        if(rind <= (int)((ir[i] - rind)/dim))
            triunz++;
    }
    blockptr->numentries = triunz;
    
    //Allocate A Memory (+1s due to Fortran index)
    blockptr->entries   = (double*)malloc((triunz+1)*sizeof(double));
    if(blockptr->entries == NULL) {
        sprintf(msgbuf,"Error allocating entries memory for A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    blockptr->iindices  = (unsigned*)malloc((triunz+1)*sizeof(unsigned));
    if(blockptr->iindices == NULL) {
        sprintf(msgbuf,"Error allocating iindices memory for A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    blockptr->jindices  = (unsigned*)malloc((triunz+1)*sizeof(unsigned));
    if(blockptr->jindices == NULL) {
        sprintf(msgbuf,"Error allocating jindices memory for A[%d, block %d]",con_num,block);
        mexErrMsgTxt(msgbuf);
    }
    
    #ifdef DEBUG
        mexPrintf("A[%d, block %d] dim: %d, triunz: %d: \n",con_num,block,dim,triunz);
    #endif 
        
    //Fill in A values, building row and column indices as we go
    ind = 1; //Fortran start
    for(i=jc[con_num];i<jc[con_num+1];i++) {
        //Row & Col Index
        rind = ir[i] % dim;
        cind = (int)((ir[i] - rind)/dim);
        //If row <= col, in upper tri
        if(rind <= cind) {
            //Fill In Indices & Val (+1 Fortran)
            blockptr->iindices[ind] = rind+1;
            blockptr->jindices[ind] = cind+1;
            blockptr->entries[ind] = -pr[i];
            
            #ifdef DEBUG
                mexPrintf("  iind[%d] = %d, jind[%d] = %d, entr[%d] = %f\n",ind,rind+1,ind,cind+1,ind,pr[i]);
            #endif
            ind++;
        }
    }      
}

//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    size_t ndec, i;
    
    //Correct number of inputs
    if(nrhs < 3)
        mexErrMsgTxt("You must supply at least 3 arguments to csdp (f, A, b)"); 
    
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
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && !mxIsStruct(pOPTS))
        mexErrMsgTxt("The options argument must be a structure");
    
    //Get Sizes
    ndec = mxGetNumberOfElements(pF);
    
    //Check linear constraints
    if(!mxIsEmpty(pA)) {
        //Check types
        if(!mxIsSparse(pA) || mxIsComplex(pA) || !mxIsDouble(pA))
            mexErrMsgTxt("Linear Constraints A must be real, sparse, double matrix");
        if(mxIsSparse(pB) || mxIsComplex(pB) || !mxIsDouble(pB))
            mexErrMsgTxt("Linear Constraints b must be real, dense, double column vector");
        //Check sizes
        if(mxGetN(pA) != ndec)
            mexErrMsgTxt("Linear Constraints A does not have the same number of columns as decision variables!");
        if(mxGetM(pA) != mxGetNumberOfElements(pB))
            mexErrMsgTxt("Linear Constraints A and b do not have the same number of rows!");        
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
}

//Check cone for correct dims and data type
void checkCone(const mxArray *cone, size_t ndec, size_t block)
{
    char msgbuf[1024];
    double Msq;
    //Check type
    if(!mxIsSparse(cone) || mxIsComplex(cone) || !mxIsDouble(cone)) {
        sprintf(msgbuf,"Cone %d must be a real, sparse, double matrix",block);
        mexErrMsgTxt(msgbuf);
    }
    //Check we can make a square matrix from M
    Msq = sqrt((double)mxGetM(cone));
    if(floor(Msq) != Msq) {
        sprintf(msgbuf,"Cone %d does not contain a number of rows that can be converted to a square matrix! CSDPMEX expects all elements.",block);
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
int user_exit(int n, int k, struct blockmatrix C, double *a, double dobj, double pobj, double constant_offset, struct constraintmatrix *con, struct blockmatrix X, double *y, struct blockmatrix Z, struct paramstruc params)
{   
    double evaltime;
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
  
    //Check max time
    if(evaltime > maxtime)
        return CSDP_MAX_TIME;

    //Check ctrl-c
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting CSDP...\n\n");
        return CSDP_USER_TERMINATION; //terminate
    }

    if(printLevel>1 && citer) {
        //Display heading if % 20 iters
        if(citer==1 || citer % 20 == 0)
            mexPrintf("Iter   Time[s]    PP Objective      DD Objective\n");
            
        //Display parameters
        mexPrintf("%-3d   %6.2f   %16.8e  %16.8e\n",citer,evaltime,pobj,dobj);
        mexEvalString("drawnow;"); //flush draw buffer
    }
    citer++;
    //Return ok
    return(0);
}

//Print Solver Information
void printSolverInfo()
{       
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" CSDP: A C Library for Semidefinite Programming [v%s]\n",CSDP_VERSION);
    PRINT_BUILD_INFO;
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/Csdp/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    #endif

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}