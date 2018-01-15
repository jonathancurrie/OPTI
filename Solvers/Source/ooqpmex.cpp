/* OOQPMEX - A MATLAB MEX Interface to OOQP
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2018
 * www.inverseproblem.co.nz
 */

/* Based in parts on ooqp_mex.c supplied with OOQP */

#include "mex.h"
#include "string.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "MehrotraSolver.h"
#include "GondzioSolver.h"
#include "QpGenSparsePardiso.h"
#include "QpGenSparseMa27.h"
#include "QpGenSparseMa57.h"
#include "OoqpMonitor.h"
#include "OoqpVersion.h"
#include "Status.h"
#include <exception>
#include <time.h>
#include <algorithm>
#include "opti_util.h"

#ifdef LINK_MKL
    #include "mkl.h"
#endif

//Linear Solvers
#define USE_PARDISO 0
#define USE_MA57 1
#define USE_MA27 2
//QP Algorithms
#define MEHROTRA 0
#define GONDZIO 1
//Extra Return Args
#define MAX_TIME_EXCEEDED -1
#define USER_EXITED -2

//Ctrl-C Detection
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);

//Define for debugging
//#define DEBUG

//Argument Enums (in expected order of arguments)
enum {eH, eF, eA, eRL, eRU, eAEQ, eBEQ, eLB, eUB, eOPTS};                   
//PRHS Defines    
#define pH      prhs[eH]
#define pF      prhs[eF]
#define pA      prhs[eA]
#define pRL     prhs[eRL]
#define pRU     prhs[eRU]
#define pAEQ    prhs[eAEQ]
#define pBEQ    prhs[eBEQ]
#define pLB     prhs[eLB]
#define pUB     prhs[eUB]
#define pOPTS   prhs[eOPTS]

using namespace std;

//Function Prototypes
void printSolverInfo();
void checkInputs(const mxArray *prhs[], int nrhs);
void lower(char *str);
void sparseTranspose(mwIndex *sJc, mwIndex *sIr, double *sPr, int *dJc, int *dIr, double *dPr, mwIndex nnz, mwIndex sM, mwIndex sN);

//Message Class
class mexPrinter : public OoqpMonitor {
    public:
        virtual void doIt( Solver * solver, Data * data, Variables * vars,
					 Residuals * resids,
					 double alpha, double sigma,
					 int i, double mu, 
                     int status_code,
					 int level );
};
//Print Handler
void mexPrinter::doIt( Solver * solver, Data * data, Variables * vars,
					 Residuals * resids,
					 double alpha, double sigma,
					 int i, double mu, 
                     int status_code,
					 int level )
{
    try
    {
        if(level < 2)
        {
            if(i == 1 || !(i%10))
                mexPrintf(" iter   duality gap           mu     resid norm\n");

            mexPrintf("%5d     %9.3g    %9.3g      %9.3g\n",i,resids->dualityGap(),mu,resids->residualNorm());
            mexEvalString("drawnow;"); //flush draw buffer
        }
    }
    catch (std::exception& error) 
    {
        mexErrMsgTxt(error.what());
    }
}
//Custom Status Class
class DerivedStatus : public Status {
    public:
        DerivedStatus(int maxIter, double maxTime, clock_t start);
        virtual int doIt( Solver * solver, Data * data,
                        Variables * vars, Residuals * resids,
                        int i, double mu, int stage );
    private:
        int _maxIter;
        double _maxTime;
        clock_t _start;
};
//Status Handler
DerivedStatus::DerivedStatus(int maxIter, double maxTime, clock_t start) : Status()
{
    _maxIter = maxIter;
    _maxTime = maxTime;
    _start = start;
}
int DerivedStatus::doIt( Solver * solver, Data * data,
                         Variables * vars, Residuals * resids,
                         int i, double mu, int stage )
{
    //Check for Ctrl-C
    if (utIsInterruptPending()) {
        utSetInterruptPending(false); /* clear Ctrl-C status */
        mexPrintf("\nCtrl-C Detected. Exiting OOQP...\n\n");
        return USER_EXITED;
    }
    //Check for MaxIts
    if(i>=_maxIter)
        return (int)MAX_ITS_EXCEEDED;    
    //Check for MaxTime
    if(((double)(clock()-_start))/CLOCKS_PER_SEC >= _maxTime)
        return MAX_TIME_EXCEEDED;
    //Check Normal Termination Tests    
    return solver->defaultStatus(data,vars,resids,min(i,5),mu,stage); //note fixed iter as we check above
}   

//Main Function
void mexFunction( int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[] ) 
{
    //Input Args
    double *H = NULL, *f = NULL, *A = NULL, *rl = NULL, *ru = NULL, *Aeq = NULL, *beq = NULL, *lb = NULL, *ub = NULL;
    
    //Return Args
    double *x, *fval, *exitflag, *iter;
    
    //Options
    int maxIter = 100;
    int printLevel = 0;    
    int linSolver = USE_MA57;
    int algorithm = GONDZIO;
    double objbias = 0;
    double maxTime = 1000;
    
    //Sparse Indicing
    double *Ht, *At, *Aeqt;
    int *iH_ir, *iH_jc, *iA_ir, *iA_jc, *iAeq_ir, *iAeq_jc;
    
    //Problem Size
    size_t ndec, neq, nin;
    mwSize nnzH, nnzA, nnzAeq;
    
    //Infinite Check Indexing
    char *ilb = NULL, *iub = NULL, *irl = NULL, *iru = NULL;    
    //Local Copies
    double *qp_rl = NULL, *qp_ru = NULL, *qp_lb  = NULL , *qp_ub = NULL;
    
    //Internal Vars
    int i, err;
    const char *fnames[4] = {"pi","y","phi","gamma"};
    char msgbuf[1024];
    mxArray *d_pi, *d_y, *d_phi, *d_gam;
    
    //Check # Inputs
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else {
            sprintf(msgbuf,"%d.%02d.%02d",OOQPVERSIONMAJOR,OOQPVERSIONMINOR,OOQPVERSIONPATCHLEVEL);
            plhs[0] = mxCreateString(msgbuf);
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }
    //Thorough Check
    checkInputs(prhs,nrhs); 
    
    //Get pointers to Input variables
    f = mxGetPr(pF);
    if(!mxIsEmpty(pA)) {
        rl = mxGetPr(pRL);
        ru = mxGetPr(pRU);
    }
    if(nrhs > eAEQ && !mxIsEmpty(pAEQ)) {
        beq = mxGetPr(pBEQ);
    }
    if(nrhs > eLB && !mxIsEmpty(pLB))
        lb = mxGetPr(pLB);
    if(nrhs > eUB && !mxIsEmpty(pUB))
        ub = mxGetPr(pUB);
    
    //Get options if specified
    if(nrhs > eOPTS) {
        if(mxGetField(pOPTS,0,"objbias"))
            objbias = *mxGetPr(mxGetField(pOPTS,0,"objbias"));
        if(mxGetField(pOPTS,0,"display"))
            printLevel = (int)*mxGetPr(mxGetField(pOPTS,0,"display"));
        if(mxGetField(pOPTS,0,"maxiter"))
            maxIter = (int)*mxGetPr(mxGetField(pOPTS,0,"maxiter"));
        if(mxGetField(pOPTS,0,"maxtime"))
            maxTime = *mxGetPr(mxGetField(pOPTS,0,"maxtime"));
        if(mxGetField(pOPTS,0,"linear_solver")) {
            if(mxIsChar(mxGetField(pOPTS,0,"linear_solver"))) {
                char *str = mxArrayToString(mxGetField(pOPTS,0,"linear_solver"));
                lower(str);
                if(!strcmp(str,"pardiso"))
                    linSolver = USE_PARDISO;
                else if(!strcmp(str,"ma57"))
                    linSolver = USE_MA57;
                else if(!strcmp(str,"ma27"))
                    linSolver = USE_MA27;
                else
                    mexErrMsgTxt("Unknown linear solver - options are 'pardiso', 'ma57' or 'ma27'");
                mxFree(str);
            }
            else
                linSolver = (int)*mxGetPr(mxGetField(pOPTS,0,"linear_solver"));
        }            
        if(mxGetField(pOPTS,0,"algorithm")) {
            if(mxIsChar(mxGetField(pOPTS,0,"algorithm"))) {
                char *str = mxArrayToString(mxGetField(pOPTS,0,"algorithm"));
                lower(str);
                if(!strcmp(str,"gondzio"))
                    algorithm = GONDZIO;
                else if(!strcmp(str,"mehrotra"))
                    algorithm = MEHROTRA;
                else
                    mexErrMsgTxt("Unknown algorithm - options are 'gondzio' or 'mehrotra'");
                mxFree(str);
            }
            else
                algorithm = (int)*mxGetPr(mxGetField(pOPTS,0,"algorithm"));            
        }
        CheckOptiVersion(pOPTS);
    }
    
    //Check Linear Solver is Available
    switch(linSolver)
    {
        case USE_PARDISO:
            #if !defined(LINK_PARDISO) && !defined(LINK_MKL_PARDISO)
                mexErrMsgTxt("PARDISO is selected as the linear solver but is not available in this build");
            #endif
            break;
        case USE_MA57:
            #if !defined(LINK_ML_MA57) && !defined(LINK_MA57)
                mexErrMsgTxt("MA57 is selected as the linear solver but is not available in this build");
            #endif
            break;
        case USE_MA27:
            #ifndef LINK_MA27
                mexErrMsgTxt("MA27 is selected as the linear solver but is not available in this build");
            #endif
            break;
        default:
            mexErrMsgTxt("Unknown linear solver selected");
    }

    //Get sizes
    ndec = mxGetNumberOfElements(pF); //f    
    nin = mxGetM(pA);           //A
    if(nrhs > eAEQ && !mxIsEmpty(pAEQ))
        neq = mxGetM(pAEQ);     //Aeq
    else
        neq = 0;
    if(mxIsEmpty(pH))
        nnzH = 0;
    else
        nnzH = mxGetJc(pH)[mxGetN(pH)];
    if(mxIsEmpty(pA))
        nnzA = 0;
    else
        nnzA = mxGetJc(pA)[mxGetN(pA)];
    if(nrhs <= eAEQ || mxIsEmpty(pAEQ))
        nnzAeq = 0;
    else
        nnzAeq = mxGetJc(pAEQ)[mxGetN(pAEQ)];
    
    //Create Outputs
    plhs[0] = mxCreateDoubleMatrix(ndec,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nin+neq,1,mxREAL);
    x = mxGetPr(plhs[0]); 
    fval = mxGetPr(plhs[1]); 
    exitflag = mxGetPr(plhs[2]);    
    iter = mxGetPr(plhs[3]);
    //Optional Outputs
    if(nlhs > 4) {    
        plhs[4] = mxCreateStructMatrix(1,1,4,fnames);
        d_pi = mxCreateDoubleMatrix(nin,1, mxREAL);
        d_y = mxCreateDoubleMatrix(neq,1, mxREAL);        
        d_phi = mxCreateDoubleMatrix(ndec,1, mxREAL);
        d_gam = mxCreateDoubleMatrix(ndec,1, mxREAL);
    }
    try
    {
        //QP Variables
        QpGenSparseMa27 *qp27 = NULL;
        QpGenSparseMa57 *qp57 = NULL;
        QpGenSparsePardiso *qpPD = NULL;
        QpGenData *prob = NULL;
        QpGenVars *vars = NULL;
        QpGenResiduals *resid = NULL;
        MehrotraSolver *sm = NULL;
        GondzioSolver *sg = NULL;
        mexPrinter *printer; //deleted automatically (I think)
        
        //Bounds Local Vectors
        ilb = (char*) mxCalloc(ndec, sizeof(char));
        iub = (char*) mxCalloc(ndec, sizeof(char));
        qp_lb = (double*)mxCalloc(ndec, sizeof(double));        
        qp_ub = (double*)mxCalloc(ndec, sizeof(double));
        //Copy lb if exists
        if(nrhs > eLB && !mxIsEmpty(pLB)) {
            for(i=0;i<ndec;i++) {
                //Create Finite lb Index Vectors + Copy Finite Values
                if(mxIsFinite(lb[i])) {
                    ilb[i] = 1;
                    qp_lb[i] = lb[i];
                }
                else {
                    ilb[i] = 0;
                    qp_lb[i] = 0.0;
                }
            }
        }
        //Else fill lb with 0s
        else {
            for(i=0;i<ndec;i++) {
                ilb[i] = 0;
                qp_lb[i] = 0.0;
            }            
        }
        //Copy ub if exists
        if(nrhs > eUB && !mxIsEmpty(pUB)) {    
            for(i=0;i<ndec;i++) {
                //Create Finite ub Index Vectors + Copy Finite Values
                if(mxIsFinite(ub[i])) {
                    iub[i] = 1;
                    qp_ub[i] = ub[i];
                }
                else {
                    iub[i] = 0;
                    qp_ub[i] = 0.0;
                }
            }
        }
        //Else fill ub with 0s
        else {
            for(i=0;i<ndec;i++) {
                iub[i] = 0;
                qp_ub[i] = 0.0;
            }            
        }
        
        //Copy Linear Inequalities RHS if exist
        if(nin > 0) {
            irl = (char*) mxCalloc(nin, sizeof(char));
            iru = (char*) mxCalloc(nin, sizeof(char));
            qp_rl = (double*)mxCalloc(nin, sizeof(double));
            qp_ru = (double*)mxCalloc(nin, sizeof(double));
            for( i = 0; i < nin; i++ ) {
                //Create Finite rl Index Vectors + Copy Finite Values
                if(mxIsFinite(rl[i])) {
                    irl[i] = 1;
                    qp_rl[i] = rl[i];
                }
                else {
                    irl[i] = 0;
                    qp_rl[i] = 0.0;
                }
                //Create Finite ru Index Vectors + Copy Finite Values
                if(mxIsFinite(ru[i])) {
                    iru[i] = 1;
                    qp_ru[i] = ru[i];
                }
                else {
                    iru[i] = 0;
                    qp_ru[i] = 0.0;
                }
            }
        }
        
        //QP H Matrix
        if(!mxIsEmpty(pH)) {
            iH_jc = (int*)mxCalloc(mxGetM(pH)+1,sizeof(int));
            iH_ir = (int*)mxCalloc(nnzH,sizeof(int)); 
            Ht = (double*)mxCalloc(nnzH,sizeof(double));
            //Transpose Matrix & convert to int32
            sparseTranspose(mxGetJc(pH), mxGetIr(pH), mxGetPr(pH), iH_jc, iH_ir, Ht, nnzH, mxGetM(pH), mxGetN(pH));           
        }
        else {
            iH_ir = (int*)mxCalloc(1,sizeof(int)); 
            iH_jc = (int*)mxCalloc(ndec+1,sizeof(int)); iH_jc[ndec] = 0;
            Ht = (double*)mxCalloc(1,sizeof(int));
        }
        
        //Linear Inequality Constraints A Matrix
        if(!mxIsEmpty(pA)) {
            //Allocate transpose memory            
            iA_jc = (int*)mxCalloc(mxGetM(pA)+1,sizeof(int));
            iA_ir = (int*)mxCalloc(nnzA,sizeof(int));             
            At = (double*)mxCalloc(nnzA,sizeof(double));
            //Transpose Matrix & convert to int32
            sparseTranspose(mxGetJc(pA), mxGetIr(pA), mxGetPr(pA), iA_jc, iA_ir, At, nnzA, mxGetM(pA), mxGetN(pA)); 
        }
        else {
            iA_ir = (int*)mxCalloc(1,sizeof(int)); 
            iA_jc = (int*)mxCalloc(1,sizeof(int)); iA_jc[0] = 0;
            At = (double*)mxCalloc(1,sizeof(double));
        }
        
        //Linear Equality Constraints Aeq Matrix
        if(nrhs > eAEQ && !mxIsEmpty(pAEQ)) {
            //Allocate transpose memory     
            iAeq_jc = (int*)mxCalloc(mxGetM(pAEQ)+1,sizeof(int));
            iAeq_ir = (int*)mxCalloc(nnzAeq,sizeof(int));              
            Aeqt = (double*)mxCalloc(nnzAeq,sizeof(double));
            //Transpose Matrix & convert to int32
            sparseTranspose(mxGetJc(pAEQ), mxGetIr(pAEQ), mxGetPr(pAEQ), iAeq_jc, iAeq_ir, Aeqt, nnzAeq, mxGetM(pAEQ), mxGetN(pAEQ)); 
        }
        else {
            iAeq_ir = (int*)mxCalloc(1,sizeof(int)); 
            iAeq_jc = (int*)mxCalloc(1,sizeof(int)); iAeq_jc[0] = 0;
            Aeqt = (double*)mxCalloc(1,sizeof(double));            
        }
        
        #ifdef DEBUG
            mexPrintf("\n\nProblem properties\n");
            for(i = 0; i < ndec; i++) {
                mexPrintf("%d: f %f lb %f ilb %d ub %f iub %d\n",i,f[i],qp_lb[i],ilb[i],qp_ub[i],iub[i]);
            }
            mexPrintf("\n");
            for(i = 0; i < nin; i++) {
                mexPrintf("%d: rl %f irl %d ru %f iru %d\n",i,qp_rl[i],irl[i],qp_ru[i],iru[i]);
            }
            mexPrintf("\n");
            for(i = 0; i < neq; i++) {
                mexPrintf("%d: beq %f\n",i,beq[i]);
            }
            mexPrintf("\n");
            H = mxGetPr(pH);
            for(i = 0; i < nnzH; i++) {
                if(i <= mxGetM(pH))
                    mexPrintf("%d: Hjc %3d Hir %3d Hpr %f\n",i,iH_jc[i], iH_ir[i], H[i]);
                else
                    mexPrintf("%d:         Hir %3d Hpr %f\n",i,iH_ir[i], H[i]);
            }
            mexPrintf("\n");
            for(i = 0; i < nnzA; i++) {
                if(i <= mxGetM(pA))
                    mexPrintf("%d: Ajc %3d Air %3d Apr %f\n",i,iA_jc[i], iA_ir[i], At[i]);
                else
                    mexPrintf("%d:         Air %3d Apr %f\n",i,iA_ir[i], At[i]);
            }
            mexPrintf("\n");
            for(i = 0; i < nnzAeq; i++) {
                if(i <= mxGetM(pAEQ))
                    mexPrintf("%d: Aeqjc %3d Aeqir %3d Aeqpr %f\n",i,iAeq_jc[i], iAeq_ir[i], Aeqt[i]);
                else
                    mexPrintf("%d:           Aeqir %3d Aeqpr %f\n",i,iAeq_ir[i], Aeqt[i]);
            }
        #endif
            
        //Create Problem, fill in data and make variables
        switch(linSolver)
        {
            #if defined(LINK_PARDISO) || defined(LINK_MKL_PARDISO)
            case USE_PARDISO:
                qpPD = new QpGenSparsePardiso((int)ndec,(int)neq,(int)nin,(int)nnzH,(int)nnzAeq,(int)nnzA);
                prob = (QpGenData*) qpPD->makeData( f, iH_jc, iH_ir, Ht, qp_lb, ilb, qp_ub, iub, iAeq_jc, iAeq_ir, Aeqt, beq, iA_jc, iA_ir, At, qp_rl, irl, qp_ru, iru);
                vars  = (QpGenVars*)qpPD->makeVariables(prob);
                resid = (QpGenResiduals*)qpPD->makeResiduals(prob);
                break;
            #endif
                
            #ifdef LINK_MA27
            case USE_MA27:
                qp27 = new QpGenSparseMa27((int)ndec,(int)neq,(int)nin,(int)nnzH,(int)nnzAeq,(int)nnzA);
                prob = (QpGenData*) qp27->makeData( f, iH_jc, iH_ir, Ht, qp_lb, ilb, qp_ub, iub, iAeq_jc, iAeq_ir, Aeqt, beq, iA_jc, iA_ir, At, qp_rl, irl, qp_ru, iru);
                vars  = (QpGenVars*)qp27->makeVariables(prob);
                resid = (QpGenResiduals*)qp27->makeResiduals(prob);
                break;
            #endif
            
            #if defined(LINK_MA57) || defined(LINK_ML_MA57)
            case USE_MA57:
                qp57 = new QpGenSparseMa57((int)ndec,(int)neq,(int)nin,(int)nnzH,(int)nnzAeq,(int)nnzA);
                prob = (QpGenData*) qp57->makeData( f, iH_jc, iH_ir, Ht, qp_lb, ilb, qp_ub, iub, iAeq_jc, iAeq_ir, Aeqt, beq, iA_jc, iA_ir, At, qp_rl, irl, qp_ru, iru);
                vars  = (QpGenVars*)qp57->makeVariables(prob);
                resid = (QpGenResiduals*)qp57->makeResiduals(prob);
                break;
            #endif
        }        
        //Make Solver
        switch(algorithm)
        {
            case MEHROTRA:
                switch(linSolver)
                {
                    case USE_PARDISO:
                        sm = new MehrotraSolver(qpPD,prob);
                        break;
                    case USE_MA27:
                        sm = new MehrotraSolver(qp27,prob);
                        break;
                    case USE_MA57:
                        sm = new MehrotraSolver(qp57,prob);
                        break;
                }                
                break;
            case GONDZIO:
                switch(linSolver)
                {
                    case USE_PARDISO:
                        sg = new GondzioSolver(qpPD,prob);
                        break;
                    case USE_MA27:
                        sg = new GondzioSolver(qp27,prob);
                        break;
                    case USE_MA57:
                        sg = new GondzioSolver(qp57,prob);
                        break;
                }  
                break;
            default:
                throw exception("Unknown algorithm!");
        }
        //Assign Status Handler (need more info from solver class really..)
        DerivedStatus *ctrlCStatus = new DerivedStatus(maxIter, maxTime, clock());
        switch(algorithm)
        {
            case MEHROTRA:
                sm->useStatus(ctrlCStatus);
                break;
            case GONDZIO:
                sg->useStatus(ctrlCStatus);
                break;
        }
        
        //Setup Options
        if(printLevel > 0) {
            if(printLevel > 1) { //iter print
                printer = new mexPrinter();
                switch(algorithm)
                {
                    case MEHROTRA:
                        sm->addMonitor( printer );
                        break;
                    case GONDZIO:
                        sg->addMonitor( printer );
                        break;
                }                
            }
            
            //Print Header
            char verStr[1024], algStr[128], linStr[128];
            getOoqpVersionString( verStr, 1024);
            switch(algorithm)
            {
                case MEHROTRA: sprintf(algStr,"Mehrotra"); break;
                case GONDZIO: sprintf(algStr,"Gondzio"); break;
            }
            switch(linSolver)
            {
                case USE_PARDISO: sprintf(linStr,"PARDSIO"); break;
                case USE_MA27: sprintf(linStr,"MA27"); break;
                case USE_MA57: sprintf(linStr,"MA57"); break;
            }
            mexPrintf("\n------------------------------------------------------------------\n");
            mexPrintf(" This is %s\n Authors: E. Michael Gertz, Stephen J. Wright\n\n",verStr);            
            mexPrintf(" This run uses the %s Solver with %s\n\n Problem Properties:\n",algStr,linStr);
            mexPrintf(" # Decision Variables:     %6d [%d nz]\n # Inequality Constraints: %6d [%d nz]\n # Equality Constraints:   %6d [%d nz]\n",ndec,nnzH,nin,nnzA,neq,nnzAeq);
            
            if(printLevel > 1)
                mexPrintf("------------------------------------------------------------------\n");
            mexEvalString("drawnow;");
        }
        
        //Solve QP
        try
        {
            switch(algorithm)
            {
                case MEHROTRA:
                    err = sm->solve(prob,vars,resid);
                    break;
                case GONDZIO:
                    err = sg->solve(prob,vars,resid);
                    break;
            }             
        }
        catch(...)
        {
            mexWarnMsgTxt("Error solving problem with OOQP");
            return;
        }

        //Assign variables
        *exitflag = err;       
        vars->x->copyIntoArray(x);
        *fval = prob->objectiveValue(vars)+objbias;
        switch(algorithm)
        {
            case MEHROTRA:
                *iter = (double)sm->iter;
                break;
            case GONDZIO:
                *iter = (double)sg->iter;
                break;
        }        
        //Assign Dual Solution
        if(nlhs > 4) {            
            vars->pi->copyIntoArray(mxGetPr(d_pi));     //dual row in
            mxSetField(plhs[4],0,fnames[0],d_pi);         
            vars->y->copyIntoArray(mxGetPr(d_y));       //dual row eq
            mxSetField(plhs[4],0,fnames[1],d_y); 
            vars->phi->copyIntoArray(mxGetPr(d_phi));   //dual col upper
            mxSetField(plhs[4],0,fnames[2],d_phi);
            vars->gamma->copyIntoArray(mxGetPr(d_gam)); //dual col lower
            mxSetField(plhs[4],0,fnames[3],d_gam);
        }

        if(printLevel > 0){
            //Termination Detected
            switch(err)
            {
                case SUCCESSFUL_TERMINATION: mexPrintf("\n *** SUCCESSFUL TERMINATION ***\n"); break;
                case MAX_ITS_EXCEEDED: mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); break;
                case MAX_TIME_EXCEEDED: mexPrintf("\n *** MAXIMUM TIME REACHED ***\n"); break;
                case USER_EXITED: mexPrintf("\n *** USER EXITED ***\n"); break;
                case INFEASIBLE: mexPrintf("\n *** TERMINATION: PROBABLY INFEASIBLE ***\n"); break;
                default: mexPrintf("\n *** TERMINATION: STATUS UNKNOWN ***\n"); break;                   
            }
            if(err == SUCCESSFUL_TERMINATION || err == MAX_ITS_EXCEEDED || err == MAX_TIME_EXCEEDED)
                mexPrintf(" Final Objective Value: %9.3g\n In %3.0f iterations\n",*fval,*iter);
            
            mexPrintf("------------------------------------------------------------------\n\n");
        }
        
         /* Free any scratch arrays */
        if( iru != NULL ) mxFree( iru ); iru = NULL;
        if( irl != NULL ) mxFree( irl ); irl = NULL;
        if( iub != NULL ) mxFree( iub ); iub = NULL;
        if( ilb != NULL ) mxFree( ilb ); ilb = NULL;   
        
        //Free Local Copies
        if( qp_ru != NULL ) mxFree( qp_ru ); qp_ru = NULL;
        if( qp_rl != NULL ) mxFree( qp_rl ); qp_rl = NULL;
        if( qp_ub != NULL ) mxFree( qp_ub ); qp_ub = NULL;
        if( qp_lb != NULL ) mxFree( qp_lb ); qp_lb = NULL;
                        
        //Free sparse memory
        mxFree(iH_ir);  mxFree(iH_jc); mxFree(Ht);
        mxFree(iA_ir);  mxFree(iA_jc); mxFree(At);
        mxFree(iAeq_ir); mxFree(iAeq_jc); mxFree(Aeqt);
        
        //Free up classes
        if( qp27 != NULL) delete qp27; qp27 = NULL;
        if( qp57 != NULL) delete qp57; qp57 = NULL;
        if( qpPD != NULL) delete qpPD; qpPD = NULL;
        if( prob != NULL) delete prob; prob = NULL;
        if( vars != NULL) delete vars; vars = NULL;
        if( resid != NULL) delete resid; resid = NULL;
        if( sm != NULL) delete sm; sm = NULL;
        if( sg != NULL) delete sg; sg = NULL;
    }
    catch(exception& e) //Unfortunately still crashes matlab...
    {
        char msgbuf[1024];
        sprintf(msgbuf,"Caught OOQP Error: %s",e.what());
        mexErrMsgTxt(msgbuf);   
    }
    catch(...)
    {
        mexErrMsgTxt("Fatal OOQP Error");
    }
}


//Check all inputs for size and type errors
void checkInputs(const mxArray *prhs[], int nrhs)
{
    int i;
    double *beq;
    size_t ndec, nin, neq;
    
    //Correct number of inputs
    if(nrhs < 5)
        mexErrMsgTxt("You must supply at least 5 arguments to OOQP (H, f, A, rl, ru)");  

    //Check we have an objective
    if(mxIsEmpty(pH) && mxIsEmpty(pF))
        mexErrMsgTxt("You must supply an objective function!");
    //Ensure we have f
    if(mxIsEmpty(pF))
        mexErrMsgTxt("You must supply f for this interface");
    
    //Check options is a structure
    if(nrhs > eOPTS && !mxIsEmpty(pOPTS) && !mxIsStruct(pOPTS))
        mexErrMsgTxt("The options argument must be a structure!");
    
    //Get Sizes
    ndec = mxGetNumberOfElements(pF);
    nin = mxGetM(pA); 
    if(nrhs > eAEQ) {
        neq = mxGetM(pAEQ); 
        if(nrhs <= eBEQ)
            mexErrMsgTxt("You must supply both Aeq and beq");
    }
    else {
        neq = 0;
    }
    
    //Check Constraint Pairs
    if(nin && mxIsEmpty(pRL))
        mexErrMsgTxt("rl is empty!");
    if(nin && mxIsEmpty(pRU))
        mexErrMsgTxt("ru is empty!");
    
    //Check Sparsity (only supported in A, AEQ and H)
    if(mxIsSparse(pF) || mxIsSparse(pRL) || mxIsSparse(pRU))
        mexErrMsgTxt("f, rl, and ru must be dense");
    if(nrhs > eBEQ && mxIsSparse(pBEQ))
        mexErrMsgTxt("beq must be dense");
    if(!mxIsEmpty(pH) && !mxIsSparse(pH))
        mexErrMsgTxt("H must be a sparse matrix");
    if(!mxIsEmpty(pA) && !mxIsSparse(pA))
        mexErrMsgTxt("A must be a sparse matrix");
    if(nrhs > eAEQ && !mxIsEmpty(pAEQ) && !mxIsSparse(pAEQ))
        mexErrMsgTxt("Aeq must be a sparse matrix");
    
    //Check Sizes
    if(!mxIsEmpty(pH) && ((mxGetM(pH) != ndec) || (mxGetN(pH) != ndec)))
        mexErrMsgTxt("H has incompatible dimensions");
    if(nin) {
        if(mxGetN(pA) != ndec)
            mexErrMsgTxt("A has incompatible dimensions");
        if(mxGetNumberOfElements(pRL) != nin)
            mexErrMsgTxt("rl has incompatible dimensions");
        if(mxGetNumberOfElements(pRU) != nin)
            mexErrMsgTxt("ru has incompatible dimensions");
    }
    if(neq) {
        if(mxGetN(pAEQ) != ndec)
            mexErrMsgTxt("Aeq has incompatible dimensions");
        if(mxGetNumberOfElements(pBEQ) != neq)
            mexErrMsgTxt("beq has incompatible dimensions");
        //Check for infinite equalities
        beq = mxGetPr(pBEQ);
        for(i=0;i<neq;i++) {
            if(mxIsInf(beq[i]) || mxIsNaN(beq[i]))
                mexErrMsgTxt("beq cannot contain Inf or NaN!");
        }
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
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" OOQP: Object Orientated Quadratic Programming [v%d.%02d.%02d]\n",OOQPVERSIONMAJOR,OOQPVERSIONMINOR,OOQPVERSIONPATCHLEVEL);     
    PRINT_BUILD_INFO;
    mexPrintf("  - (C) 2001 University of Chicago\n");
    mexPrintf("  - Source available from: http://pages.cs.wisc.edu/~swright/ooqp/\n\n");
    
    #if defined(LINK_MKL) || defined(LINK_MA27) || defined(LINK_MA57) || defined(LINK_MKL_PARDISO)
        mexPrintf(" This binary is statically linked to the following software:\n");
    #endif
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    #endif
    #ifdef LINK_MKL_PARDISO
        mexPrintf("  - Intel MKL PARDISO [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);
    #endif
    #ifdef LINK_MA27
        mexPrintf("  - HSL MA27 (This Binary MUST NOT BE REDISTRIBUTED)\n");
    #endif        
    #ifdef LINK_MA57
        mexPrintf("  - HSL MA57 (This Binary MUST NOT BE REDISTRIBUTED)\n");
        #ifdef LINK_METIS
            mexPrintf("  - MeTiS [v4.0.3] Copyright University of Minnesota\n");
        #endif
    #endif
    #ifdef LINK_ML_MA57
        mexPrintf("\n This binary is dynamically linked to the following software:\n");
        mexPrintf("  - MA57   [v3.0] (Included as part of the MATLAB distribution)\n");
    #endif

    mexPrintf("\n MEX Interface J.Currie 2013 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
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

//Sparse Matrix Transpose code
void sparseTranspose(mwIndex *sJc, mwIndex *sIr, double *sPr, int *dJc, int *dIr, double *dPr, mwIndex nnz, mwIndex sM, mwIndex sN)
{
    mwIndex i, j, ind = 0, index;    
            
    //Transpose Dimensions
    mwIndex dN = sM;
    mwIndex dM = sN;
    
    //Create working memory
    int *inc = (int*)mxCalloc(sM,sizeof(int));
    int *lmem = (int*)mxCalloc(nnz,sizeof(int));
    
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

	for(i = 0; i < sN; i++)  //build full missing triple
		for(j = 0; j < (sJc[i+1]-sJc[i]); j++)
			lmem[ind++] = (int)i;
		
	for(i = 0; i < nnz; i++) {
		ind = sIr[i];
		index = (dJc[ind]) + inc[ind]++;
		dIr[index] = lmem[i];    //new 'Ir' from generated sparse triple	
        dPr[index] = sPr[i];
 	}
    //Free local memory
    mxFree(inc);
    mxFree(lmem);
}