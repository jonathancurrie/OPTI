/* GSLMEX_NLS - A MATLAB MEX Interface to GSL NLS Solver
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "gslmex.h"
#include "opti_util.h"

#include <gsl_blas.h>
#include <gsl_multifit_nlinear.h>

#include <time.h>

// Inputs/Outputs
#define pPROB   prhs[0]
#define pXSOL   plhs[0]
#define pFEVAL  plhs[1]
#define pEFLAG  plhs[2]
#define pSTATS  plhs[3]

// Extra return codes
#define GSL_USER_EXIT   (-5)
#define GSL_MAX_TIME    (-6)
#define GSL_MAX_FEVAL   (-7)
#define GSL_NO_SOL      (-8)

// MATLAB data
#define FLEN 128 /* max length of user function name */
typedef struct matlab_data_s
{
    char *f = "feval";
    mxArray *plhs[1];
    mxArray *prhs[2];
    mxArray *prhs_g[2];
    size_t ndec     = 0;
    size_t ndata    = 0;
    double *ydata   = nullptr;
    int nfeval      = 0;
} matlab_data_t;

typedef struct iterfun_data_s
{
    char *f = "feval";
    mxArray *plhs[1];
    mxArray *prhs[4];
    size_t ndec   = 0;
    int verbose   = 0;
    bool enabled  = false;
} iterfun_data_t;

// Local Function Prototypes
int fun(const gsl_vector* x, void* data, gsl_vector* f);
int grad(const gsl_vector* x, void* data, gsl_matrix* J);
void iterFun(const size_t iter, void* data, const gsl_multifit_nlinear_workspace *gslWs);

// Solve a NLS using GSL
void gslSolveNLS(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs)
{
    // Check the input arguments
    
    // Determine sizes
    size_t ndec  = Mex::getFieldNumElem(pPROB, "x0"); mxGetNumberOfElements(mxGetField(pPROB,0,"x0"));
    size_t ndata = Mex::getFieldNumElem(pPROB, "ydata"); mxGetNumberOfElements(mxGetField(pPROB,0,"ydata"));
    
     //Create Outputs
    pXSOL  = Mex:createDoubleMatrix(ndec, 1); mxCreateDoubleMatrix(ndec,1, mxREAL);
    pFEVAL = Mex:createDoubleMatrix(1, 1); mxCreateDoubleMatrix(1,1, mxREAL);
    pEFLAG = Mex:createDoubleMatrix(1, 1); mxCreateDoubleMatrix(1,1, mxREAL);
    
    pSTATS = mxCreateStructMatrix(1,1, 5, fieldNames);
    
    
    
    plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(ndec, ndec, mxREAL);
    double* x = mxGetPr(pXSOL); 
    double* fval = mxGetPr(pFEVAL); 
    double* exitflag = mxGetPr(pEFLAG);    
    double* iter = mxGetPr(plhs[3]);
    double* nfeval = mxGetPr(plhs[4]);
    double* ngeval = mxGetPr(plhs[5]);
    double* covar = mxGetPr(plhs[6]);
    
    // Create the callback memory
    matlab_data_t mlData;
    mlData.prhs[0] = mxGetField(prhs[0],0,"fun");
    mlData.prhs[1] = mxCreateDoubleMatrix(ndec, 1, mxREAL); //x0
    mlData.ndec    = ndec;
    mlData.ndata   = ndata;
    mlData.ydata   = mxGetPr(mxGetField(prhs[0],0,"ydata"));
    
    // Assign parameters
    gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
    
    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
    
    // Get options if specified
    double xtol = 1e-7;
    double gtol = 1e-7;
    double ftol = 1e-7;
    int maxiter = 1000;
    int verbose = 0;
    int maxfeval = 10000;
    double maxtime = 1000;
    mxArray* opts = mxGetField(prhs[0],0,"options");
    if (opts != nullptr && !mxIsEmpty(opts))
    {
        getDoubleOption(opts, "tolafun", &ftol);
        getIntegerOption(opts, "maxiter", &maxiter);
        getIntegerOption(opts, "display", &verbose);
        getIntegerOption(opts, "maxfeval", &maxfeval);
        getDoubleOption(opts, "maxtime", &maxtime);
    }
    
    // Assign optional iteration callback
    iterfun_data_t iterData;
    iterData.verbose = verbose;
    iterData.ndec    = ndec;
    
    // Create solver memory
    gsl_multifit_nlinear_workspace* gslWs = gsl_multifit_nlinear_alloc(T, &params, ndata, ndec);
    gsl_vector* x_gsl = gsl_vector_alloc(ndec);
    
    // Copy in x0
    double* x0 = mxGetPr(mxGetField(prhs[0],0,"x0"));
    for (size_t i = 0; i < ndec; i++)
    {
        gsl_vector_set(x_gsl, i, x0[i]);
    }
    
    // Create the function definition structure
    gsl_multifit_nlinear_fdf fdf;
    fdf.f = fun;
    fdf.df = nullptr;
    fdf.fvv = nullptr;
    fdf.n = ndata;
    fdf.p = ndec;
    fdf.params = &mlData;    
    
    // Initialize the solver
    
    int status = gsl_multifit_nlinear_init (x_gsl, &fdf, gslWs);
    if (status != GSL_SUCCESS)
    {
        gsl_multifit_nlinear_free(gslWs);
        gsl_vector_free(x_gsl);
        MEX_ERR("Error initializing GSL solver! Error code: %d\n", status);
        return;
    }
    
    //Print Header
    if(verbose > 0) 
    {
        mexPrintf("\n------------------------------------------------------------------\n");        
        mexPrintf(" This is GSL v%s\n",GSL_VERSION);           
        mexPrintf(" Algorithm: %s using %s\n", gsl_multifit_nlinear_name(gslWs), gsl_multifit_nlinear_trs_name(gslWs));
        mexPrintf(" Problem Properties:\n");
        mexPrintf(" # Decision Variables:     %4d\n",ndec);
        mexPrintf(" # Data Points:            %4d\n",ndata);
        mexPrintf("------------------------------------------------------------------\n");
    }
    
    // Call GSL Solver
    status  = 0;
    int info    = 0;
    size_t niter = 0;
    clock_t startTime = clock();
    do
    {
        // Solver Iteration
        status = gsl_multifit_nlinear_iterate (gslWs);
        if (status == GSL_ENOPROG && iter == 0)
        {
              status = GSL_NO_SOL;
              break;
        }

        // Call the callback function 
        if ((iterData.verbose > 1) || (iterData.enabled == true))
        {
            iterFun(++niter, &iterData, gslWs);
        }
        
        // Test for convergence
        status = gsl_multifit_nlinear_test(xtol, gtol, ftol, &info, gslWs);
        
        // Check for Ctrl-C
        if (utIsInterruptPending()) 
        {
            utSetInterruptPending(false); /* clear Ctrl-C status */
            mexPrintf("\nCtrl-C Detected. Exiting GSL...\n\n");
            status = GSL_USER_EXIT;
            break;
        }

        // Check for max fevals
        if (mlData.nfeval > maxfeval)
        {
            status = GSL_MAX_FEVAL;
            break;
        }

        //Get Execution Time
        clock_t now = clock();
        double evaltime = ((double)(now - startTime))/CLOCKS_PER_SEC;
        // Check for max time
        if(evaltime > maxtime)
        {
            status = GSL_MAX_TIME;
            break;
        }
    }
    while ((status == GSL_CONTINUE) && (niter < maxiter));
    
    // Check return status
    if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
        info = status;
        status = GSL_SUCCESS;
    }
    if ((niter >= maxiter) && (status != GSL_SUCCESS))
    {
        status = GSL_EMAXITER;
    }
    
    // Assign solution
    gsl_vector* x_sol = gsl_multifit_nlinear_position(gslWs);
    memcpy(x, x_sol->data, ndec * sizeof(double));     
    *exitflag = static_cast<double>(status);
    *iter     = static_cast<double>(gsl_multifit_nlinear_niter(gslWs));
    *nfeval   = static_cast<double>(fdf.nevalf);
    *ngeval   = static_cast<double>(fdf.nevaldf);
    // Assign fval
    gsl_vector* f = gsl_multifit_nlinear_residual(gslWs);
    gsl_blas_ddot(f, f, fval);
    
    if (status == GSL_SUCCESS)
    {
        // Assign covariance
        gsl_matrix* J = gsl_multifit_nlinear_jac(gslWs);
        gsl_matrix* covar_gsl = gsl_matrix_alloc (ndec, ndec);
        gsl_multifit_nlinear_covar (J, 0.0, covar_gsl);
        // Copy to ML memory
        memcpy(covar, covar_gsl->data, ndec*ndec*sizeof(double));
        gsl_matrix_free(covar_gsl);
    }
    
    //Print Footer
    if(verbose > 0)
    { 
        switch(status)
        {
            case GSL_SUCCESS:
            {
                if (info == 0)
                {
                    mexPrintf("\n *** SUCCESSFUL TERMINATION: Gradient Convergence ***\n"); 
                }
                else
                {
                    mexPrintf("\n *** SUCCESSFUL TERMINATION: Step Size Convergence ***\n"); 
                }               
                mexPrintf(" Final SSE: %12.5g\n In %3.0f iterations\n",*fval,*iter);
                break;
            }
            case GSL_EMAXITER:
            {
                mexPrintf("\n *** MAXIMUM ITERATIONS REACHED ***\n"); 
                break;
            }
            case GSL_ENOPROG:
            {
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: No Further Progress Possible ***\n");
                break;
            }
            case GSL_MAX_TIME:
            {
                mexPrintf("\n *** MAXIMUM SOLVING TIME REACHED ***\n"); 
                break;
            }
            case GSL_MAX_FEVAL:
            {
                mexPrintf("\n *** MAXIMUM FUNCTION EVALUATIONS REACHED ***\n"); 
                break;
            }
            case GSL_NO_SOL:
            {
                mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: No Progress in First Iteration ***\n");
                break;
            }
            default:
            {
                mexPrintf("\n *** TERMINATION: ROUTINE ERROR (Code %d) ***\n",(int)*exitflag); 
                break;
            }
        }
        mexPrintf("------------------------------------------------------------------\n\n");
    }
    
    // Free solver memory
    gsl_multifit_nlinear_free(gslWs);
    gsl_vector_free(x_gsl);
}


// Objective function
static int fun(const gsl_vector* x, void* data, gsl_vector* f)
{
    // Get the ML data
    matlab_data_t* mlData = static_cast<matlab_data_t*>(data);
    
    // Copy in x
    mlData->plhs[0] = nullptr;
    memcpy(mxGetPr(mlData->prhs[1]), x->data, mlData->ndec * sizeof(double));

    // Call MATLAB
    if(mexCallMATLAB(1, mlData->plhs, 2, mlData->prhs, mlData->f) != 0)
    {
        mexErrMsgTxt("Error calling Objective Function!");
    }
    mlData->nfeval++;
    
    // Copy Objective
    double *fval = mxGetPr(mlData->plhs[0]);
    for (size_t i = 0; i < mlData->ndata; i++)
    {
        gsl_vector_set(f, i, fval[i] - mlData->ydata[i]);
    }
    
    // Clean up
    mxDestroyArray(mlData->plhs[0]);
    
    return GSL_SUCCESS;
}

// Iteration Callback Function
static void iterFun(const size_t iter, void* data, const gsl_multifit_nlinear_workspace *gslWs)
{
    // Get iteration callback data
    iterfun_data_t* iterData = static_cast<iterfun_data_t*>(data);
    
    // Calculate norm
    gsl_vector *f = gsl_multifit_nlinear_residual(gslWs);
    double norm = gsl_blas_dnrm2(f);
    double sse  = norm*norm;

    // Optional print output
    if (iterData->verbose > 1)
    {
        // Computer some statistics
        double rcond  = 0.0;
        gsl_multifit_nlinear_rcond(&rcond, gslWs);
        
        // Display header every 10 iters
        if(!(iter % 10))
        {
            mexPrintf(" iter          cond(J)                sse\n");
        }    
        mexPrintf("%5d         %8.5g       %12.5g\n", iter, 1.0 / rcond, sse);
        mexEvalString("drawnow;"); //flush draw buffer  
    }
    
    //Iteration Callback
    if(iterData->enabled)
    {
        int niter = static_cast<int>(iter);
        gsl_vector* x = gsl_multifit_nlinear_position(gslWs);
        
        // Assign iteration data
        iterData->plhs[0] = NULL;        
        memcpy(mxGetData(iterData->prhs[1]), &niter, sizeof(int));
        memcpy(mxGetPr(iterData->prhs[2]), &sse, sizeof(double));
        memcpy(mxGetPr(iterData->prhs[3]), x->data, iterData->ndec * sizeof(double));
        
        // Call MATLAB
        if (mexCallMATLAB(1, iterData->plhs, 4, iterData->prhs, iterData->f) != 0)
        {
            mexErrMsgTxt("Error calling Callback Function!");
        }

        //Collect return argument
        bool stop = *(bool*)mxGetData(iterData->plhs[0]);
        
        // Clean up
        mxDestroyArray(iterData->plhs[0]);
        
        //Warn user stop not implemented
        if(stop)
        {
            mexWarnMsgTxt("GSL NLS does not implement the stop feature of iterfun");
        }
    }
}
    
    