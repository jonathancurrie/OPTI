/* FILTERSDMEX - A MATLAB MEX Interface to FILTERSD
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.inverseproblem.co.nz
 */

/* This MEX file contains two interfaces; one for the dense version of
 * filterSD and one for the sparse version. The default is a dense
 * interface. To enable the sparse interface, pass the preprocessor
 * SPARSEVER.                   
 */

/* This code also contains a snippet from Peter Carbonetto's IPOPT
 * interface, released under the EPL 
 */

#include "mex.h"
#include <time.h>
#include <string.h>
#include "opti_util.h"

#ifdef LINK_MKL
    #include "mkl.h"
#endif

//FILTERSD Version
#define FILTERSD_VERSION "1.0"

//Enable for debug print out
//#define DEBUG
//Enable for debug print out on sparse transpose
//#define DEBUGTR

//Function handle structure
#define FLEN 128 /* max length of user function name */
#define MAXRHS 2 /* max nrhs for user function */
typedef struct {
     char f[FLEN], c[FLEN];
     mxArray *f_prhs[MAXRHS], *f_plhs[1];
     mxArray *c_prhs[MAXRHS], *c_plhs[1];
     int f_xrhs, f_nrhs;
     int c_xrhs, c_nrhs;
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

//Argument Enums (in expected order of arguments)
#ifdef SPARSEVER
    enum {eFUN, eGRAD, eX0, eLB, eUB, eNLCON, eNLJAC, eNLJACSTR, eCL, eCU, eOPTS};
    #define pNLJACSTR  prhs[eNLJACSTR]
#else
    enum {eFUN, eGRAD, eX0, eLB, eUB, eNLCON, eNLJAC, eCL, eCU, eOPTS};        
#endif            
//PRHS Defines    
#define pFUN    prhs[eFUN]
#define pGRAD   prhs[eGRAD]
#define pX0     prhs[eX0]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pNLCON  prhs[eNLCON]
#define pNLJAC  prhs[eNLJAC]
#define pCL     prhs[eCL]
#define pCU     prhs[eCU]
#define pOPTS   prhs[eOPTS]
//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
void sparseTranspose(mwIndex *sJc, mwIndex *sIr, double *sPr, mwIndex *dJc, mwIndex *dIr, double *dPr, mwIndex nnz, mwIndex sM, mwIndex sN);
// FILTERSD routine implemented in Fortran.
extern void FILTERSD (int *n, int *m, double *x, double *al, double *f, double *fmin, char *cstype,
                      double *bl, double *bu, double *ws, int *lws, double *v, int *nv, int *maxa, 
                      int *maxla, int *maxu, int *maxiu, int *kmax, int *maxg, double *rho,
                      double *htol, double *rgtol, int *maxit, int *iprint, int *nout, int *ifail);    

//Global Variables (not passing via user_data)
user_function_data fun;
user_function_data grad;
iter_fun_data iterF;
int printLevel, citer, maxfev;
//Max Time data
double maxtime;
clock_t start, end;
//Jacobian Structure Sparse Indexing Variables (needed to compare to returned Jacs)
#ifdef SPARSEVER
    mwIndex *str_jc, *str_ir, str_nz;
    mwSize str_m, str_n;
    mwIndex *lmem = NULL, *inc = NULL; //transpose variables
    mxArray *jacT = NULL;
    mwIndex *jacT_jc, *jacT_ir;
    double *jacT_pr;
#endif

//Common blocks in FORTRAN
extern struct {
    int kk, ll, kkk, lll;
    int mxws_, mxlws_;
} WSC;

extern struct {
    double ainfty, ubd;
    int mlp, mxf;
} DEFAULTC;

extern struct {
    double fxd,alc;
    int m_,iph,last1,next1,nx,nx1,nal,nal1,naal,naal1,nxd,nxd1,ncx,ncx1,ncxd,ncxd1,nla1;
} FUNCTC;
    
extern struct {
    double dnorm, h, hJt, hJ;
    int ipeq, k, itn, nft, ngt;
} STATSC;

#ifdef SPARSEVER
    extern struct {
        int ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,neb,neb1,nprof;
        int lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1;
        int nq,nq1,nr,nr1,ny,ny1,nz,nz1,lv,lv1,le,le1;
    } SCHURC;
    extern struct {
        int mc, mxmc;
    } REFACTORC;
#else
    extern struct {
        int ns,ns1,nt,nt1,nu,nu1,mx1,lc,lc1,li,li1;
    } DENSEC;
    extern struct {
        int mxm1;
    } MXM1C;    
#endif

// Function definitions. 
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
    //Expected calling form
    //filtersd(fun,grad,x0,lb,ub,nlcon,nljac,nljacstr,cl,cu,opts)
    
    //Outputs Args
    double *x, *fval, *exitflag, *nfval, *ngval, *niter, *al;
    const char *fnames[3] = {"niter","nfval","ngval"};
    
    //Internal Vars
    int i, j;
    double *lb, *llb, *ub, *lub, *x0, *cl, *cu;
    char msgbuf[1024];
    double evaltime;
    #ifdef SPARSEVER
        mwIndex *jc, *ir, k;
        mxArray *aT;
    #endif

    //FILTERSD Vars
    int ndec = 0;
    int ncon = 0;
    int maxg = 7;           //max reduced gradient vectors to store
    int maxu = 0;           //length of user workspace
    int maxiu = 0;          //length of user integer workspace
    int maxa = 0;           //max entries in Jacobian (set by gradients) [modified below]
    int maxla = 0;          //locations required for sparse matrix indices and pointers 
    int kmax;
    int ifail = 0;          //failure code
    int maxit = 1000;       //maximum iterations
    int iprint = 0;         //print value
    double rho = 1e1;       //initial trust region radius
    double htol = 1e-6;     //tolerance for sum of constraint infeasibilities
    double rgtol = 1e-4;    //tolerance in reduced gradient L2 norm
    double fmin;
    double *xint;           //internal x
    char *cstype;           //character worksapce
    double v[7];            //stores nv Ritz values
    int nv = 1;             //set for no previous info
    //Workspace Variables  
    int mxmc,mxm1,np,nm,kk,kkk,ll;
    int wsmem = 100;//10000000;
    int lwsmem = 50;//5000000;
    double *ws;
    int *lws;
        
    //Check for no inputs
    if (nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString("1.0");          
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }

    //Check user inputs
    checkInputs(prhs,nrhs);
    
    //Set Defaults
    citer = 1;                  //feval counter
    iterF.enabled = false;      //not iteration callback
    printLevel = 0;             //no message output
    maxfev = 10000;             //maximum function evaluations
    maxit = 1000;               //maximum iterations
    maxtime = 1000;             //maximum time

    //Get Sizes
    ndec = (int)mxGetNumberOfElements(pX0);
    if(nrhs > eNLCON)
        ncon = (int)mxGetNumberOfElements(pCL);
        
    //Set some defaults
    maxg = min(ndec,6)+1;       //max gradients to store
    kmax = ndec;                //must be <= ndec (max dimension of null space)
    
    //Attempt to calculate required WS and LWS memory 
    #ifdef SPARSEVER
        if(nrhs > eNLCON && !mxIsEmpty(pNLCON)) {
            jc = mxGetJc(pNLJACSTR);               
            str_nz = jc[mxGetN(pNLJACSTR)];
        }
        else
            str_nz = 0;       
        maxa = ((int)str_nz)+ndec;  //max entries in sparse Jacobian   
        maxla = maxa+ncon+3;        //max entries in filterSD sparse format
    #else
        maxa = ndec*(ncon+1);       //max entries in Jacobian (note +1 includes gradient)                         
        maxla = 1;
    #endif
    //Calc common intermediate totals            
    mxmc = min(ndec,500);       //filterSD default  
    mxm1 = min(ncon+1,ndec);    //non-trivial block of basis matrix 
    nm = ndec+ncon;
    np = 2*mxmc + 5*ndec + (mxmc*(mxmc + 1))/2 + mxmc*mxmc + 2;      
    kk = 6*ncon + 2*maxa + maxu + DEFAULTC.mlp + 2*DEFAULTC.mxf + 8*ndec;                
    kkk = nm+ndec+maxg*(maxg+1)/2+maxg*(kmax+5);    
    ll = maxla + maxiu + DEFAULTC.mlp + nm;
    //Determine workspace requirements based on sparsity
    #ifdef SPARSEVER
        //WS Workspace
        wsmem = np + kk + kkk;        
        //LWS Workspace
        lwsmem = ll + WSC.lll + 2*mxmc + 8*ndec + nm + 2;
        //Remember this is below min, so increase by 50% (opti heuristic)
        wsmem = (int)(1.5*(double)wsmem); 
    #else
        //WS Workspace
        wsmem = kk + kkk + mxm1*(mxm1+1)/2 + 3*ndec + mxm1;
        //LWS Workspace
        lwsmem = ll + WSC.lll + ndec + mxm1 + nm;
    #endif
    #ifdef DEBUG
        mexPrintf("kk: %d, kkk: %d, ll: %d, np: %d, mxmc: %d, mxm1 %d\n",kk,kkk,ll,np,mxmc,mxm1);
        mexPrintf("Calculated ws Workspace: %d, lws Workspace: %d\n\n",wsmem,lwsmem);
    #endif
    
    //Set FilterSD Workspace Memory based on Above Calc
    WSC.mxws_  = wsmem;
    WSC.mxlws_ = lwsmem;
    //Allocate memory for workspace
    ws = mxCalloc(WSC.mxws_,sizeof(double));
    lws = mxCalloc(WSC.mxlws_,sizeof(int));
        
    //Get Objective Function Handle
    if (mxIsChar(pFUN)) {
        CHECK(mxGetString(pFUN, fun.f, FLEN) == 0,"error reading objective name string");
        fun.f_nrhs = 1;
        fun.f_xrhs = 0;
    } else {
        fun.f_prhs[0] = (mxArray*)pFUN;
        strcpy(fun.f, "feval");
        fun.f_nrhs = 2;
        fun.f_xrhs = 1;
    }
    fun.f_prhs[fun.f_xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    //Get Gradient Function Handle
    if (mxIsChar(pGRAD)) {
        CHECK(mxGetString(pGRAD, grad.f, FLEN) == 0,"error reading gradient name string");
        grad.f_nrhs = 1;
        grad.f_xrhs = 0;
    } else {
        grad.f_prhs[0] = (mxArray*)pGRAD;
        strcpy(grad.f, "feval");
        grad.f_nrhs = 2;
        grad.f_xrhs = 1;
    }   
    grad.f_prhs[grad.f_xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    //Get x0
    x0 = mxGetPr(pX0);
    
    //Setup Lower Bounds  (note also includes CL) 
    lb = mxCalloc(ndec+ncon,sizeof(double));
    if(nrhs > eLB && !mxIsEmpty(pLB)) {
        llb = mxGetPr(pLB);
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(llb[i])) 
                lb[i] = llb[i];
            else
                lb[i] = -DEFAULTC.ainfty;
        }
    }
    else {
        for(i=0;i<ndec;i++)
            lb[i] = -DEFAULTC.ainfty;
    }
    //Setup Upper Bounds (note also includes CU)
    ub = mxCalloc(ndec+ncon,sizeof(double));
    if(nrhs > eUB && !mxIsEmpty(pUB)) {
        lub = mxGetPr(pUB);
        for(i=0;i<ndec;i++) {
            if(!mxIsInf(lub[i])) 
                ub[i] = lub[i];
            else
                ub[i] = DEFAULTC.ainfty;
        }
    }
    else {
        for(i=0;i<ndec;i++)
            ub[i] = DEFAULTC.ainfty;
    }
    
    //Get nlcon if specified
    if(nrhs > eNLCON && !mxIsEmpty(pNLCON)) {
        //Get Nonlinear Constraint Function Handle
        if (mxIsChar(pNLCON)) {
            CHECK(mxGetString(pNLCON, fun.c, FLEN) == 0,"error reading nonlinear constraint name string");
            fun.c_nrhs = 1;
            fun.c_xrhs = 0;
        } else {
            fun.c_prhs[0] = (mxArray*)pNLCON;
            strcpy(fun.c, "feval");
            fun.c_nrhs = 2;
            fun.c_xrhs = 1;
        }
        fun.c_prhs[fun.c_xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
        //Get Jacobian Function Handle
        if (mxIsChar(pNLJAC)) {
            CHECK(mxGetString(pNLJAC, grad.c, FLEN) == 0,"error reading jacobian name string");
            grad.c_nrhs = 1;
            grad.c_xrhs = 0;
        } else {
            grad.c_prhs[0] = (mxArray*)pNLJAC;
            strcpy(grad.c, "feval");
            grad.c_nrhs = 2;
            grad.c_xrhs = 1;
        }   
        grad.c_prhs[grad.c_xrhs] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0 
        //Get Lower Row Bounds
        cl = mxGetPr(pCL);
        for(i=0,j=ndec;i<ncon;i++,j++) {
            if(!mxIsInf(cl[i])) 
                lb[j] = cl[i];
            else
                lb[j] = -DEFAULTC.ainfty;
        }
        //Get Upper Row Bounds
        cu = mxGetPr(pCU);
        for(i=0,j=ndec;i<ncon;i++,j++) {
            if(!mxIsInf(cu[i])) 
                ub[j] = cu[i];
            else
                ub[j] = DEFAULTC.ainfty;
        }
        
        //Process Jacobian Structure if sparse
        #ifdef SPARSEVER
            //These variables are BEFORE transpose
            str_m = mxGetM(pNLJACSTR);      
            str_n = mxGetN(pNLJACSTR);      
            jc = mxGetJc(pNLJACSTR);        
            ir = mxGetIr(pNLJACSTR);        
            str_nz = jc[str_n]; 
            //Create transpose working memory
            inc = mxCalloc(str_m,sizeof(mwIndex));
            lmem = mxCalloc(str_nz,sizeof(mwIndex)); 
            //Create a new copy of the Jacobian Structure and transpose it for filterSD
            aT = mxCreateSparse(str_n,str_m,str_nz,mxREAL);
            str_jc = mxGetJc(aT); str_ir = mxGetIr(aT);
            sparseTranspose(jc, ir, mxGetPr(pNLJACSTR), str_jc, str_ir, mxGetPr(aT), str_nz, str_m, str_n);            
            //Generate lws vector (pointer locations for sparse entries - but not really pointers)
            //Remember FORTRAN indexing +1 (although starts at zero anyway)
            i = 0;
            lws[i++] = ((int)str_nz)+ndec+1;  //total entries + 1
            //Gradient Start
            for(j=1;j<=ndec;j++)
                lws[i++] = j;               //gradient always dense
            //Jacobian Start
            for(j=0;j<(int)str_nz;j++)
                lws[i++] = ((int)str_ir[j])+1;//row starts (col once transposed)
            lws[i++] = 1;                   //start of gradient
            lws[i++] = ndec+1;              //start of jacobian
            k = ndec+1;
            for(j=0;j<((int)str_m)-1;j++) {  
                k += str_jc[j+1] - str_jc[j];   //column (row once transposed) groups
                lws[i++] = (int)k;
            }
            lws[i] = ((int)str_nz)+ndec+1;        //as with initial value
            //Create a sparse matrix to use for transposing Jacobian into in gradient callback
            jacT = mxCreateSparse(str_n,str_m,str_nz,mxREAL);
            jacT_jc = mxGetJc(jacT); jacT_ir = mxGetIr(jacT); jacT_pr = mxGetPr(jacT);
            
            #ifdef DEBUG
                mexPrintf("nnz: %d, str_m: %d, str_n: %d, ndec: %d\n",str_nz,str_m,str_n,ndec);
                for(i=0;i<(int)str_nz+ndec+ncon+3;i++)
                    mexPrintf("lws[%d] = %d\n",i,lws[i]);
            #endif
        #endif
    
    }
    else {
        fun.c_nrhs = 0; //indicate not present
        grad.c_nrhs = 0; 
        #ifdef SPARSEVER
            //Fill in gradient starts    
            i = 0;
            lws[i++] = +ndec+1;             //total entries + 1
            //Gradient Start
            for(j=1;j<=ndec;j++)
                lws[i++] = j;               //gradient always dense
            lws[i++] = 1;                   //start of gradient
            lws[i++] = ndec+1;              //start of jacobian
        #endif
    }
    
    //Get Options if specified
    if(nrhs > eOPTS) {
        if(mxGetField(pOPTS,0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display"));
        if(mxGetField(pOPTS,0,"maxfeval"))
            maxfev = (int)*mxGetPr(mxGetField(pOPTS,0,"maxfeval"));
        if(mxGetField(pOPTS,0,"maxiter"))
            maxit = (int)*mxGetPr(mxGetField(pOPTS,0,"maxiter"));
        if(mxGetField(pOPTS,0,"maxtime"))
            maxtime = *mxGetPr(mxGetField(pOPTS,0,"maxtime"));
        if(mxGetField(pOPTS,0,"rho"))
            rho = *mxGetPr(mxGetField(pOPTS,0,"rho"));
        if(mxGetField(pOPTS,0,"htol"))
            htol = *mxGetPr(mxGetField(pOPTS,0,"htol"));
        if(mxGetField(pOPTS,0,"rgtol"))
            rgtol = *mxGetPr(mxGetField(pOPTS,0,"rgtol"));
        if(mxGetField(pOPTS,0,"iterfun") && !mxIsEmpty(mxGetField(pOPTS,0,"iterfun")))
        {
            iterF.prhs[0] = (mxArray*)mxGetField(pOPTS,0,"iterfun");
            strcpy(iterF.f, "feval");
            iterF.enabled = true;  
            iterF.prhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            iterF.prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
            iterF.prhs[3] = mxCreateDoubleMatrix(ndec,1,mxREAL);
        }
        CheckOptiVersion(pOPTS);
    } 
        
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    //Statistic Structure Output
    plhs[3] = mxCreateStructMatrix(1,1,3,fnames);
    mxSetField(plhs[3],0,fnames[0],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[1],mxCreateDoubleMatrix(1,1, mxREAL));
    mxSetField(plhs[3],0,fnames[2],mxCreateDoubleMatrix(1,1, mxREAL));
    niter  = mxGetPr(mxGetField(plhs[3],0,fnames[0]));
    nfval = mxGetPr(mxGetField(plhs[3],0,fnames[1]));
    ngval   = mxGetPr(mxGetField(plhs[3],0,fnames[2]));
    //Lambda
    plhs[4] = mxCreateDoubleMatrix(ndec+ncon,1, mxREAL);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    al = mxGetPr(plhs[4]);
    
    //Allocate memory
    xint = mxCalloc(ndec+ncon,sizeof(double));
    cstype = mxCalloc(ncon,sizeof(char));
    //Initialize
    v[0] = 1;                   //no previous ritz info
    nv = 1;                     //as above    
    #ifdef SPARSEVER
        maxa = ((int)str_nz)+ndec;  //max entries in sparse Jacobian
        maxla = maxa+ncon+3;        //max entries including column starts, begin, tail and middle
    #else
        maxa = ndec*(ncon+1);       //max entries in Jacobian (note +1 includes gradient)
        maxla = 1;                  //indicates dense        
        lws[0] = ndec;              //stride
        MXM1C.mxm1 = min(ncon+1,ndec); //non-trivial block of basis matrix (dense only) [max no of general cons in active set]
    #endif
    fmin = -DEFAULTC.ainfty;    //minimum fval   
    DEFAULTC.ubd = 1e20;        //maximum initial constraint violation
    
    //Copy in x0
    memcpy(xint,x0,ndec*sizeof(double));
    
    //Print Header
    if(printLevel) {
        mexPrintf("\n------------------------------------------------------------------\n");
        #ifdef SPARSEVER
            mexPrintf(" This is FILTERSD [Sparse]\n"); 
        #else    
            mexPrintf(" This is FILTERSD [Dense]\n"); 
        #endif
        mexPrintf(" Author: Roger Fletcher\n\n");
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:        %4d\n",ndec);
        #ifdef SPARSEVER
            mexPrintf(" # Nonlinear Constraints:     %4d [%d nz]\n",ncon,str_nz);
        #else
            mexPrintf(" # Nonlinear Constraints:     %4d\n",ncon);
        #endif

        mexPrintf("------------------------------------------------------------------\n");
        mexEvalString("drawnow;");
    }

    //Start timer
    start = clock();
    
    #ifdef DEBUG 
        mexPrintf("FILTERSD INPUT ARGS: ndec: %d, ncon: %d, nv: %d, maxa: %d, maxla: %d, maxu: %d, maxiu: %d, kmax: %d, maxg: %d, maxit: %d\n",ndec,ncon,nv,maxa,maxla,maxu,maxiu,kmax,maxg,maxit);
    #endif
    
    //Call FilterSD
    FILTERSD (&ndec,&ncon,xint,al,fval,&fmin,cstype,lb,ub,ws,lws,v,&nv,&maxa,&maxla,&maxu,&maxiu,&kmax,&maxg,&rho,&htol,&rgtol,&maxit,&iprint,NULL,&ifail);       
    //Check for not enough ws or lws
    if(ifail == 7) {
        #ifdef SPARSEVER 
            if(printLevel) {
                sprintf(msgbuf,"Sparse Workspace memory was too small, ws = %d, req >> %d, lws = %d, req = %d. Reallocating...",WSC.mxws_,WSC.mxws_-SCHURC.nprof,WSC.mxlws_,WSC.ll + WSC.lll + SCHURC.le + REFACTORC.mxmc + 1);
                mexWarnMsgTxt(msgbuf);
            }         
            //Increase workspace memory to that calculated by filterSD (+ a little extra)
            WSC.mxws_ = max(WSC.mxws_,3*(WSC.mxws_-SCHURC.nprof)); //assume 3x is 'much greater than req'
            WSC.mxlws_ = max(WSC.mxlws_,WSC.ll + WSC.lll + SCHURC.le + REFACTORC.mxmc + 1)+100;
            
        #else    
            if(printLevel) {
                sprintf(msgbuf,"Dense Workspace memory was too small, ws = %d, req = %d, lws = %d, req = %d. Reallocating...",WSC.mxws_,DENSEC.ns,WSC.mxlws_,DENSEC.nt);
                mexWarnMsgTxt(msgbuf);
            }
            //Increase workspace memory to that calculated by filterSD (+ a little extra)            
            WSC.mxws_ = max(WSC.mxws_,DENSEC.ns)+100;
            WSC.mxlws_ = max(WSC.mxlws_,DENSEC.nt)+100;  
            lws[0] = ndec;                     
        #endif 
        //Reallocate MATLAB workspace memory
        mxRealloc(ws,WSC.mxws_*sizeof(double));
        mxRealloc(lws,WSC.mxlws_*sizeof(int)); //this should preserve exisitng lws
        if(ws == NULL)
            mexErrMsgTxt("Error reallocating ws (double) workspace memory");
        if(lws == NULL)
            mexErrMsgTxt("Error reallocating lws (integer) workspace memory");
        if(printLevel) {            
            sprintf(msgbuf,"Reallocated ws to %d*%d, lws to %d*%d\n\n",WSC.mxws_,sizeof(double),WSC.mxlws_,sizeof(int));
            mexWarnMsgTxt(msgbuf);
        }
        //No previous Ritz info
        v[0] = 1; nv = 1; ifail = 0;
        //Call FilterSD again
        FILTERSD (&ndec,&ncon,xint,al,fval,&fmin,cstype,lb,ub,ws,lws,v,&nv,&maxa,&maxla,&maxu,&maxiu,&kmax,&maxg,&rho,&htol,&rgtol,&maxit,&iprint,NULL,&ifail);
    }
    else if(ifail == 8) {
        //Double filter length
        DEFAULTC.mxf *= 5;
        if(printLevel)
            sprintf(msgbuf,"Filter length was too short, extending.");
        //Call FilterSD again
        FILTERSD (&ndec,&ncon,xint,al,fval,&fmin,cstype,lb,ub,ws,lws,v,&nv,&maxa,&maxla,&maxu,&maxiu,&kmax,&maxg,&rho,&htol,&rgtol,&maxit,&iprint,NULL,&ifail);
    }
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Copy results back
    memcpy(x,xint,ndec*sizeof(double));
    *exitflag = (double)ifail;
    *niter = (double)STATSC.itn;
    *nfval = (double)STATSC.nft;
    *ngval = (double)STATSC.ngt;
    
    //Print Header
    if(printLevel){            
        //Termination Detected
        switch(ifail)
        {
            //Success
            case 0:	
                mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n"); break;
            //Error
            case 2:
                mexPrintf("\n *** ERROR: Bounds on x are inconsistent ***\n"); break;
            case 5:
                mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;
            case 105:
                if(citer > maxfev) {
                    mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n");
                    *exitflag = 101;
                }
                else if(evaltime > maxtime) {
                    mexPrintf("\n *** MAXIMUM TIME REACHED ***\n");
                    *exitflag = 102;
                }
                else
                    mexPrintf("\n *** TERMINATION: USER EXITED ***\n"); break;
                    
                break;
            //Early Exit
            case 1:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Unbounded NLP: f(x) <= fmin at an htol-feasible point x ***\n"); break;
            case 3:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Local minimum of feasibility problem and h(x) > htol (nonlinear constraints are locally inconsistent) ***\n"); break;
            case 4:
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: Initial point x has h(x) > ubd (upper constraint violation limit) ***\n"); break;
            case 6:
                mexPrintf("\n *** TERMINATION: ERROR ***\n *** CAUSE: Termination with rho <= htol ***\n"); break;
            case 7:
                mexPrintf("\n *** TERMINATION: ERROR ***\n *** CAUSE: Not enough workspace memory for solver ***\n"); break;
            case 8:
                mexPrintf("\n *** TERMINATION: ERROR ***\n *** CAUSE: Insufficient space for filter ***\n"); break;
            case 9:
            case 10:
                mexPrintf("\n *** TERMINATION: ERROR ***\n *** CAUSE: Unexpected fail in LCP solver ***\n"); break;
            default:
                mexPrintf("\n *** TERMINATION: ERROR ***\n *** CAUSE: Unknown! ***\n"); break;
        }

        if(*exitflag==0)
            mexPrintf("\n Final Objective: %12.5g\n In %5d major iterations\n    %5d function evaluations\n    %5d gradient evaluations\n",*fval,STATSC.itn,STATSC.nft,STATSC.ngt);

        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    //Free memory
    mxFree(xint); mxFree(cstype);
    mxFree(lb); mxFree(ub);
    mxFree(ws); mxFree(lws);
    #ifdef SPARSEVER
        if(jacT) {mxDestroyArray(jacT); jacT=NULL;}
        if(inc) {mxFree(inc); inc=NULL;}
        if(lmem) {mxFree(lmem); lmem=NULL;}
    #endif
}

void FUNCTIONS(int *n, int *m, double *x, double *f, double *c, void *user, int *iuser, int *ncont)
{
    int stat;
    double *cval;
    bool stop = false;
    double evaltime;
    
    //Get Execution Time
    end = clock();
    evaltime = ((double)(end-start))/CLOCKS_PER_SEC;
    
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); // clear Ctrl-C status
        if(printLevel)
            mexPrintf("\nCtrl-C Detected. Exiting FILTERSD...\n\n");
        *ncont = 0; //terminate
        return;
    }       
    //Check for maxtime expiry    
    if(evaltime > maxtime) {
        if(printLevel)
            mexPrintf("\nMaximum Solver Time Exceeded. Exiting FILTERSD...\n\n");
        *ncont = 0; //terminate
        return;
    }
    //Check for maxfeval
    if(citer > maxfev) {
        if(printLevel)
            mexPrintf("\nMaximum Function Evaluations Exceeded. Exiting FILTERSD...\n\n");
        *ncont = 0; //terminate
        return;
    }
    
    //Call MATLAB for Objective
    fun.f_plhs[0] = NULL;
    memcpy(mxGetPr(fun.f_prhs[fun.f_xrhs]), x, *n * sizeof(double));    
    stat = mexCallMATLAB(1, fun.f_plhs, fun.f_nrhs, fun.f_prhs, fun.f);
    if(stat)
      mexErrMsgTxt("Error calling Objective Function!");    
    //Get Objective
    *f = *mxGetPr(fun.f_plhs[0]);;    
    // Clean up Ptr
    mxDestroyArray(fun.f_plhs[0]);

    //Call MATLAB for Nonlinear Constraints (if present)
    if(fun.c_nrhs) {
        fun.c_plhs[0] = NULL;
        memcpy(mxGetPr(fun.c_prhs[fun.c_xrhs]), x, *n * sizeof(double));    
        stat = mexCallMATLAB(1, fun.c_plhs, fun.c_nrhs, fun.c_prhs, fun.c);
        if(stat)
          mexErrMsgTxt("Error calling Nonlinear Constraint Function!");    
        //Copy in Constraint Vector
        cval = mxGetPr(fun.c_plhs[0]);
        memcpy(c,cval,*m*sizeof(double));  
        // Clean up Ptr
        mxDestroyArray(fun.c_plhs[0]);
    }
    
    //Iteration Printing
    if(printLevel > 1) {
        if(citer == 1 || !(citer%10))
            mexPrintf(" feval       time           objective\n");

        mexPrintf("%5d       %5.2f     %15.8g\n",citer,evaltime,*f);                    
        mexEvalString("drawnow;"); //flush draw buffer
    }

    //Iteration Callback
    if(iterF.enabled)
    {
        iterF.plhs[0] = NULL;
        memcpy(mxGetData(iterF.prhs[1]), &citer, sizeof(int));
        memcpy(mxGetPr(iterF.prhs[2]), f, sizeof(double));
        memcpy(mxGetPr(iterF.prhs[3]), x, *n * sizeof(double));
        stat = mexCallMATLAB(1, iterF.plhs, 4, iterF.prhs, iterF.f);
        if(stat)
            mexErrMsgTxt("Error calling Callback Function!");

        //Collect return argument
        stop = *(bool*)mxGetData(iterF.plhs[0]);
        // Clean up Ptr
        mxDestroyArray(iterF.plhs[0]);
    }
    
    //Check for iterfun terminate
    if (stop) {
        if(printLevel)
            mexPrintf("\nIterFun called Stop. Exiting FILTERSD...\n\n");
        *ncont = 0; //terminate
    }
    
    citer++;
}

void GRADIENTS(int *n, int *m, double *x, double *a, void *user, int *iuser)
{
    int stat;
    double *gval;
    #ifdef SPARSEVER
        double *sPr;
        mwIndex *sJc, *sIr, snz, i, j, c;
        bool match, matchrow, matchcol;  
    #else
        int i,j,k;
    #endif
                
    //Call MATLAB for Gradient
    grad.f_plhs[0] = NULL;
    memcpy(mxGetPr(grad.f_prhs[grad.f_xrhs]), x, *n * sizeof(double));    
    stat = mexCallMATLAB(1, grad.f_plhs, grad.f_nrhs, grad.f_prhs, grad.f);
    if(stat)
      mexErrMsgTxt("Error calling Gradient Function!");    
    //Copy in Gradient
    gval = mxGetPr(grad.f_plhs[0]);
    memcpy(a,gval,*n*sizeof(double));  
    // Clean up Ptr
    mxDestroyArray(grad.f_plhs[0]);   
    
    //Call MATLAB for Jacobian (if present)
    if(grad.c_nrhs) {
        grad.c_plhs[0] = NULL;
        memcpy(mxGetPr(grad.c_prhs[grad.c_xrhs]), x, *n * sizeof(double));    
        stat = mexCallMATLAB(1, grad.c_plhs, grad.c_nrhs, grad.c_prhs, grad.c);
        if(stat)
          mexErrMsgTxt("Error calling Jacobian Function!");   
        #ifdef SPARSEVER
            //First initialize all entries incase we don't have values for them all
            for(i=0,j=*n;i<str_nz;i++,j++)
                a[j] = 0;
            //Now extract return matrices properties
            if(!mxIsSparse(grad.c_plhs[0]))
                mexErrMsgTxt("The Jacobian returned from MATLAB is not sparse!");
            sJc = mxGetJc(grad.c_plhs[0]);
            sIr = mxGetIr(grad.c_plhs[0]);
            sPr = mxGetPr(grad.c_plhs[0]);
            snz = sJc[mxGetN(grad.c_plhs[0])];                   
            if(snz > str_nz)
                mexErrMsgTxt("The Jacobian returned from MATLAB has more non-zero entries than the passed Jacobian Structure!");
            if(mxGetN(grad.c_plhs[0]) != str_n || mxGetM(grad.c_plhs[0]) != str_m)
                mexErrMsgTxt("The size of the Jacobian returned from MATLAB was not the same as the passed Jacobian Structure size!");
            
            //Copy and transpose Jacobian into memory
            sparseTranspose(sJc, sIr, sPr, jacT_jc, jacT_ir, jacT_pr, snz, str_m, str_n);             
            //The following code is modified from Dr. Peter Carbonetto's code in sparsematrix.cpp [part of the IPOPT MEX interface] [EPL] 
            //Copy in transposed sparse matrix elements to equivalent positions in Jacobian a (after gradient)
            i = 0;  // Index of element in source.
            j = 0;  // Index of element in destination.
            for (c = 0; c < str_m; c++) {
                // Repeat for each non-zero element in the destination column.
                for ( ; j < str_jc[c+1]; j++) {
                  matchrow = (jacT_ir[i] == str_ir[j]);
                  matchcol = (i >= jacT_jc[c]) && (i < jacT_jc[c+1]);
                  match    = matchrow && matchcol;
                  a[*n+j] = match * jacT_pr[i];
                  i += match;
                }      
            }
        #else
            //Copy in Dense Jacobian, transposing as we go and tagging onto the end of the gradient
            gval = mxGetPr(grad.c_plhs[0]); k = *n;
            for(i=0;i<*m;i++)
                for(j=0;j<*n;j++)
                    a[k++] = gval[i+j**m];
        #endif
        // Clean up Ptr
        mxDestroyArray(grad.c_plhs[0]); 
    }
}

void checkInputs(const mxArray *prhs[], int nrhs)
{    
    size_t ndec;
    
    if(nrhs < 3)
        mexErrMsgTxt("You must supply at least 3 arguments to filtersd!\n\nfiltersd(fun,grad,x0)\n");
       
    //Check Types
    if(!mxIsFunctionHandle(pFUN) && !mxIsChar(pFUN))
        mexErrMsgTxt("fun must be a function handle or function name!");
    if(!mxIsEmpty(pGRAD) && (!mxIsFunctionHandle(pGRAD) && !mxIsChar(pGRAD)))
        mexErrMsgTxt("grad must be a function handle or function name!");
    if(!mxIsDouble(pX0) || mxIsComplex(pX0) || mxIsEmpty(pX0))
        mexErrMsgTxt("x0 must be a real double column vector!");
    
    //Get ndec
    ndec = mxGetNumberOfElements(pX0);
    
    //Check Bounds
    if(nrhs > eLB) {
        if(!mxIsDouble(pLB) || mxIsComplex(pLB))
            mexErrMsgTxt("lb must be a real double column vector!");
        if(nrhs > eUB && (!mxIsDouble(pUB) || mxIsComplex(pUB)))
            mexErrMsgTxt("ub must be a real double column vector!");
        //Check Sizes
        if(!mxIsEmpty(pLB) && (ndec != mxGetNumberOfElements(pLB)))
            mexErrMsgTxt("lb is not the same length as x0! Ensure they are both Column Vectors");
        if(nrhs > eUB && !mxIsEmpty(pUB) && (ndec != mxGetNumberOfElements(pUB)))
            mexErrMsgTxt("ub is not the same length as x0! Ensure they are both Column Vectors");
    }
    
    //Check Nonlinear Constraint Handle
    if(nrhs > eNLCON && !mxIsEmpty(pNLCON)) {
        if(nrhs <= eCU) //ensure we have all args req
            mexErrMsgTxt("You must supply all nonlinear arguments (including Jacobian, cl, cu and structure if sparse)");        
        if(!mxIsFunctionHandle(pNLCON) && !mxIsChar(pNLCON))
            mexErrMsgTxt("nlcon must be a function handle or function name!");
        if(!mxIsFunctionHandle(pNLJAC) && !mxIsChar(pNLJAC))
            mexErrMsgTxt("nljac must be a function handle or function name!");
        #ifdef SPARSEVER
            if(!mxIsSparse(pNLJACSTR))
                mexErrMsgTxt("nljacstr must be a sparse matrix!");
            if(mxGetM(pNLJACSTR) != mxGetNumberOfElements(pCL))
                mexErrMsgTxt("nljacstr does not have the same number of rows as constraints");
            if(mxGetN(pNLJACSTR) != ndec)
                mexErrMsgTxt("nljacstr does not have the same number of columns as variables");
        #endif           
        //Check row bounds
        if(mxIsEmpty(pCL) || mxIsEmpty(pCU))
            mexErrMsgTxt("cl and cu must be specified when supplying nonlinear constraints");
        if(!mxIsDouble(pCL) || mxIsComplex(pCL))
            mexErrMsgTxt("cl must be a real double column vector!");
        if(!mxIsDouble(pCU) || mxIsComplex(pCU))
            mexErrMsgTxt("cu must be a real double column vector!");
        //Check Sizes
        if(mxGetNumberOfElements(pCL) != mxGetNumberOfElements(pCU))
            mexErrMsgTxt("cl is not the same length as cu!");
    }

    //Check Options
    if(nrhs > eOPTS) {
        if(!mxIsStruct(pOPTS))
            mexErrMsgTxt("The specified options must be a structure!");
    }
}

#ifdef SPARSEVER
//Sparse Matrix Transpose code [bit rough - if anyone knows a better way let me know!]
void sparseTranspose(mwIndex *sJc, mwIndex *sIr, double *sPr, mwIndex *dJc, mwIndex *dIr, double *dPr, mwIndex nnz, mwIndex sM, mwIndex sN)
{
    mwIndex i, j, ind = 0, index;    
            
    //Transpose Dimensions
    mwIndex dN = sM;
    mwIndex dM = sN;
    
    //Initialize working memory
    for(i=0;i<dN;i++)
        dJc[i] = 0;
    for(i=0;i<dN;i++)
        inc[i] = 0;
            
    //Transpose and fill in
	for(i = 0; i < nnz; i++)
		dJc[sIr[i]+1]++;	//sum each row index to determine no elements in each row
	for(i = 2; i <= dN; i++)
		dJc[i] += dJc[i-1]; //cumsum to determine new 'Jc'
 
    #ifdef DEBUGTR
        for(i = 0; i <= dN; i++)
            mexPrintf("sJc[%d]: %d, dJc: %d\n",i,sJc[i],dJc[i]);
    #endif
    
	for(i = 0; i < sN; i++)  //build full missing triple
		for(j = 0; j < (sJc[i+1]-sJc[i]); j++)
			lmem[ind++] = i;
		
	for(i = 0; i < nnz; i++) {
		ind = sIr[i];
		index = (dJc[ind]) + inc[ind]++;
		dIr[index] = lmem[i];    //new 'Ir' from generated sparse triple	
        dPr[index] = sPr[i];
        #ifdef DEBUGTR
            mexPrintf("%d: index %d ind %d inc %d lmem %d\n",i,index,ind,inc[ind],lmem[i]);
        #endif
 	}
}
#endif

//Print Solver Information
void printSolverInfo()
{    
    mexPrintf("\n-----------------------------------------------------------\n");
    #ifdef SPARSEVER
        mexPrintf(" FILTERSD: FilterSD [Sparse Version] Nonlinear Optimizer [v%s]\n",FILTERSD_VERSION);
    #else
        mexPrintf(" FILTERSD: FilterSD [Dense Version] Nonlinear Optimizer [v%s]\n",FILTERSD_VERSION);
    #endif  
    PRINT_BUILD_INFO;
    mexPrintf("  - Released under the Eclipse Public License: http://opensource.org/licenses/eclipse-1.0\n");
    mexPrintf("  - Source available from: https://projects.coin-or.org/filterSD\n\n");

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}