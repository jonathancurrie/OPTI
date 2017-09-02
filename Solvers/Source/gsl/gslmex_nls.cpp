/* GSLMEX_NLS - A MATLAB MEX Interface to GSL NLS Solver
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "Include/gslmex_nls.h"
#include "opti_mex_utils.h"
#include <gsl_version.h>
#include <gsl_blas.h>

namespace opti_gsl_nls
{
// Inputs/Outputs
#define pPROB   (args.prhs[0])
#define pXSOL   (args.plhs[0])
#define pFEVAL  (args.plhs[1])
#define pEFLAG  (args.plhs[2])
#define pSTATS  (args.plhs[3])

using namespace opti_mex_utils;

// Local Function Definitions
static int fun(const gsl_vector* x, void* data, gsl_vector* f);
static int grad(const gsl_vector* x, void* data, gsl_matrix* df);
static bool iterFun(const size_t iter, void* data, const gsl_multifit_nlinear_workspace *gslWs);

// Solve a NLS using GSL
void GSLMexNLS::solve(const opti_mex_utils::OptiMexArgs& args)
{
    // Check the input arguments, will throw if there's an error
    checkInputArgs(args);
    
    // Determine sizes
    size_t ndec  = OPTIMex::getFieldNumElem(pPROB, "x0"); 
    size_t ndata = OPTIMex::getFieldNumElem(pPROB, "ydata"); 
    
     // Create Outputs
    pXSOL  = OPTIMex::createDoubleMatrix(ndec, 1); 
    pFEVAL = OPTIMex::createDoubleScalar(0.0); 
    pEFLAG = OPTIMex::createDoubleScalar(0.0);
    const char* stats[] = {"niter", "nfeval", "ngeval", "covar"};
    pSTATS = OPTIMex::createStruct(stats, 4);
    
    // Collect Pointers to Outputs
    double* x           = mxGetPr(pXSOL); 
    double* fval        = mxGetPr(pFEVAL); 
    double* exitflag    = mxGetPr(pEFLAG);    
    double* niter       = OPTIMex::createFieldScalar(pSTATS, "niter", 0.0);
    double* nfeval      = OPTIMex::createFieldScalar(pSTATS, "nfeval", 0.0);
    double* ngeval      = OPTIMex::createFieldScalar(pSTATS, "ngeval", 0.0);
    double* covar       = OPTIMex::createFieldMatrix(pSTATS, "covar", ndec, ndec);
    
    // Create the callback memory
    matlab_cb_data_t mlData;
    iterfun_cb_data_t iterData;
    OPTIMex::initMatlabCallbackData(mlData, pPROB, ndec, ndata);
    
    // Assign default solver parameters
    gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
    
    const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
    
    // Set Default Options
    double xtol     = OPTI_DEFAULT_TOL;
    double gtol     = OPTI_DEFAULT_TOL;
    double ftol     = OPTI_DEFAULT_TOL;
    int verbose     = OPTI_VERBOSE_OFF;
    OPTITermSettings termSettings("GSL_NLS");
    // Update options if specified
    if (OPTIMex::isValidField(pPROB, "options"))
    {
        mxArray* opts = OPTIMex::getField(pPROB, "options");
        OPTIMex::getDoubleOption(opts, "tolafun", ftol);
        OPTIMex::getIntegerOption(opts, "maxiter", termSettings.maxIter);
        OPTIMex::getIntegerOption(opts, "display", verbose);
        OPTIMex::getIntegerOption(opts, "maxfeval", termSettings.maxFeval);
        OPTIMex::getDoubleOption(opts, "maxtime", termSettings.maxTime);

        // Optional iteration callback
        if (OPTIMex::isValidField(opts, "iterfun"))
        {
            OPTIMex::initIterFunCallbackData(iterData, OPTIMex::getField(opts, "iterfun"), ndec);
        }
    }
    
    // Create solver memory
    gsl_multifit_nlinear_workspace* gslWs = gsl_multifit_nlinear_alloc(T, &params, ndata, ndec);
    gsl_vector* x_gsl = gsl_vector_alloc(ndec);
    
    try
    {
        // Copy in x0
        double* x0 = OPTIMex::getFieldPr(pPROB,"x0");
        for (size_t i = 0; i < ndec; i++)
        {
            gsl_vector_set(x_gsl, i, x0[i]);
        }
        
        // Create the function definition structure
        gsl_multifit_nlinear_fdf fdf;
        fdf.f       = fun;
        if (mlData.haveGrad == true)
        {
            fdf.df  = grad;
        }
        fdf.fvv     = nullptr;
        fdf.n       = ndata;
        fdf.p       = ndec;
        fdf.params  = &mlData;    
        
        // Initialize the solver
        int status = gsl_multifit_nlinear_init (x_gsl, &fdf, gslWs);
        if (status != GSL_SUCCESS)
        {
            gsl_multifit_nlinear_free(gslWs);
            gsl_vector_free(x_gsl);
            OPTIMex::error("Error initializing GSL solver! Error code: %d\n", status);
            return;
        }
        
        //Print Header
        if (verbose != OPTI_VERBOSE_OFF) 
        {
            printHeader(args.solverInfo, gslWs);
        }
        
        // Call GSL Solver
        status          = 0;
        int info        = 0;
        size_t iter     = 0;
        // Create and start the timer
        OPTITimer timer;
        do
        {
            // Solver Iteration
            status = gsl_multifit_nlinear_iterate (gslWs);
            if (status == GSL_ENOPROG && iter == 0)
            {
                status = GSL_NO_SOL;
                break;
            }
            // Increment iteration count
            ++iter;

            // Print Solver Update
            if (verbose >= OPTI_VERBOSE_ITER)
            {
                printIter(gslWs, timer.execTime());
            }

            // Call the iteration callback function 
            if (iterData.enabled == true)
            {
                if (iterFun(iter, &iterData, gslWs) == true) // User requested exit
                {
                    status = GSL_USER_EXIT;
                    break;
                }
            }

            // Check for termination
            switch (OPTIMex::checkTermination(timer, iter, mlData.nfeval, 0, termSettings))
            {
                case OPTI_TERM_USER_EXIT: { status = GSL_USER_EXIT; break; }
                case OPTI_TERM_MAX_FEVAL: { status = GSL_MAX_FEVAL; break; }
                case OPTI_TERM_MAX_ITER:  { status = GSL_EMAXITER; break; }
                case OPTI_TERM_MAX_TIME:  { status = GSL_MAX_TIME; break; }
                default:                  { status = gsl_multifit_nlinear_test(xtol, gtol, ftol, &info, gslWs); break; }
            }
        }
        while (status == GSL_CONTINUE);
        
        // Check return status
        if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
        {
            info = status; // Used for Footer
            status = GSL_SUCCESS;
        }
        
        // Assign solution
        gsl_vector* x_sol = gsl_multifit_nlinear_position(gslWs);
        memcpy(x, x_sol->data, ndec * sizeof(double));   
        // Assign fval
        gsl_vector* f = gsl_multifit_nlinear_residual(gslWs);
        gsl_blas_ddot(f, f, fval);
        // Copy Statistics  
        *exitflag = static_cast<double>(status);
        *niter    = static_cast<double>(gsl_multifit_nlinear_niter(gslWs));
        *nfeval   = static_cast<double>(fdf.nevalf);
        *ngeval   = static_cast<double>(fdf.nevaldf);
        
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
        if(verbose != OPTI_VERBOSE_OFF)
        {
            printFooter(status, info, *fval, *niter, *exitflag);        
        }
    }
    catch (...)
    {
        // Free solver memory
        gsl_multifit_nlinear_free(gslWs);
        gsl_vector_free(x_gsl);

        // Rethrow exception
        throw;
    }

    // Free solver memory
    gsl_multifit_nlinear_free(gslWs);
    gsl_vector_free(x_gsl);
}


// Objective function
static int fun(const gsl_vector* x, void* data, gsl_vector* f)
{
    // Get the ML data
    matlab_cb_data_t* mlData = static_cast<matlab_cb_data_t*>(data);

    // Call MATLAB
    double *fval = OPTIMex::callMatlabObjective(mlData, x->data);

    // Copy residual
    for (size_t i = 0; i < mlData->ndata; i++)
    {
        gsl_vector_set(f, i, fval[i] - mlData->ydata[i]);
    }
    
    // Clean up
    mxDestroyArray(mlData->plhs[0]);
    return GSL_SUCCESS;
}

// Gradient function
static int grad(const gsl_vector* x, void* data, gsl_matrix* df)
{
    // Get the ML data
    matlab_cb_data_t* mlData = static_cast<matlab_cb_data_t*>(data);

    // Call MATLAB
    double *gval = OPTIMex::callMatlabGradient(mlData, x->data);

    // Copy gradient while tranposing (row major -> col major)
    OPTIMex::copyTranspose(df->data, gval, mlData->ndata, mlData->ndec);
    
    // Clean up
    mxDestroyArray(mlData->plhs[0]);
    return GSL_SUCCESS;
}

// Iteration Callback Function
static bool iterFun(const size_t iter, void* data, const gsl_multifit_nlinear_workspace *gslWs)
{
    // Get iteration callback data
    iterfun_cb_data_t* iterData = static_cast<iterfun_cb_data_t*>(data);
    
    // Calculate norm
    gsl_vector *f = gsl_multifit_nlinear_residual(gslWs);
    double norm = gsl_blas_dnrm2(f);
    double sse  = norm*norm;
    // Access decision variables
    gsl_vector* x = gsl_multifit_nlinear_position(gslWs);

    // Call MATLAB and return stop status
    return OPTIMex::callMatlabIterFun(iterData, iter, sse, x->data);
}
    

//
// Print Routines 
//
void GSLMexNLS::printIter(const gsl_multifit_nlinear_workspace* gslWs, double execTime)
{
    // Compute some statistics
    gsl_vector *f   = gsl_multifit_nlinear_residual(gslWs);
    double norm     = gsl_blas_dnrm2(f);
    double sse      = norm*norm;
    double rcond    = 0.0;
    gsl_multifit_nlinear_rcond(&rcond, gslWs);

    const char* headerString = "  sse            cond(J)\n";
    const char* formatString = "%14.5g       %12.5g\n";

    OPTIMex::printSolverIter(gslWs->niter, execTime, headerString, formatString, sse,  1.0 / rcond);
}

void GSLMexNLS::printHeader(const OptiSolverProperties& solverInfo, const gsl_multifit_nlinear_workspace* gslWs)
{
    char algBuf[1024];
    snprintf(algBuf, 1024, "%s using %s\n", gsl_multifit_nlinear_name(gslWs), gsl_multifit_nlinear_trs_name(gslWs));
    OPTIMex::printSolverHeader(solverInfo, algBuf, gslWs->x->size, gslWs->f->size);
}

void GSLMexNLS::printFooter(int status, int info, double fval, double niter, double exitflag)
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
            mexPrintf(" Final SSE: %12.5g\n In %3.0f iterations\n",fval,niter);
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
        case GSL_USER_EXIT:
        {
            mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: User Exited ***\n"); 
            break;
        }
        case GSL_NO_SOL:
        {
            mexPrintf("\n *** TERMINATION: EARLY EXIT ***\n *** CAUSE: No Progress in First Iteration ***\n");
            break;
        }
        default:
        {
            mexPrintf("\n *** TERMINATION: ROUTINE ERROR (Code %d) ***\n",(int)exitflag); 
            break;
        }
    }
    mexPrintf("------------------------------------------------------------------\n\n");
}
    
void GSLMexNLS::checkInputArgs(const OptiMexArgs& args)
{

}

 } // namespace opti_gsl_nls