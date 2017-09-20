/* MKLJAC - A MATLAB MEX Interface to Intel MKL's DJACOBI
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "mkl.h"
#include "mklJac.h"
using namespace opti_mex_utils;

// Settings
#define OPTI_DEFAULT_MKLJAC_TOL (1e-6)

// Inputs/Outputs
#define pFUN    (args.prhs[0])
#define pX      (args.prhs[1])
#define pLEN    (args.prhs[2])
#define pTOL    (args.prhs[3])
#define pDX     (args.plhs[0])
#define pSTATUS (args.plhs[1])
#define pNFEVAL (args.plhs[2])

 // Solver Defines
#define SOLVER_NAME     ("mklJac: Intel djacobi")
#define SOLVER_LICENSE  (OptiSolverLicense::MKL)
#define SOLVER_LINK     ("https://software.intel.com/en-us/mkl")

// Local Function Declarations
static void mklFun(MKL_INT *m, MKL_INT *n, double *x, double *f, void *data);

//
// Main Function
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    OptiSolverProperties solver(SOLVER_NAME, __INTEL_MKL__, __INTEL_MKL_MINOR__, __INTEL_MKL_UPDATE__, SOLVER_LINK, SOLVER_LICENSE);
    OptiMexArgs args(solver, prhs, nrhs, plhs, nlhs);

    // Check if just requesting version/solver info
    if (OPTIMex::checkNoArgs(args) == true)
    {
        return; 
    }

    // Perform the differentation & return
    MKLJac::differentiate(args);
}  


//
// Main Differentiation Method
//
void MKLJac::differentiate(const opti_mex_utils::OptiMexArgs& args)
{
    bool haveSize = false;
    bool haveTol  = false;
    // Check the input arguments
    checkInputArgs(args, haveSize, haveTol);

    // Determine size of x (n)
    double* x    = mxGetPr(pX);
    MKL_INT lenX = static_cast<MKL_INT>(OPTIMex::getNumElem(pX));

    // Create the callback memory
    matlab_cb_data_t mlData;
    OPTIMex::initMatlabCallbackFun(mlData, pFUN, lenX);

    // Determine size of f (m)
    MKL_INT lenF = 0;
    if (haveSize == true) // length supplied
    {
        lenF = static_cast<MKL_INT>(OPTIMex::getDoubleScalar(pLEN));
    }
    else  // Dummy call for length
    {
        mlData.plhs[0] = nullptr;
        mexCallMATLAB(1, mlData.plhs, 2, mlData.prhs_f, mlData.f);
        lenF = static_cast<MKL_INT>(OPTIMex::getNumElem(mlData.plhs[0]));
        mxDestroyArray(mlData.plhs[0]);
    }            

    // Get tolerance, if specified
    double tol = OPTI_DEFAULT_MKLJAC_TOL;
    if (haveTol == true)
    {
        tol = OPTIMex::getDoubleScalar(pTOL);
    }

    // Create outputs
    pDX     = OPTIMex::createDoubleMatrix(lenF, lenX);
    pSTATUS = OPTIMex::createDoubleScalar(0.0);
    pNFEVAL = OPTIMex::createDoubleScalar(0.0);
    double* dx = mxGetPr(pDX);
    double* status = mxGetPr(pSTATUS);
    double* nfeval = mxGetPr(pNFEVAL);

    // Perform the differentiation
    if (djacobix(mklFun, &lenX, &lenF, dx, x, &tol, &mlData) == TR_SUCCESS)
    {
        *status = 1.0;
    }
    else
    {
        *status = 0.0;
    }

    // Assign numFevals
    *nfeval = static_cast<double>(mlData.nfeval);
}

//
// MATLAB Callback
//
static void mklFun(MKL_INT *m, MKL_INT *n, double *x, double *f, void *data)
{
    // Get the ML data
    matlab_cb_data_t* mlData = static_cast<matlab_cb_data_t*>(data);
    // Call MATLAB
    double *fval = OPTIMex::callMatlabObjective(mlData, x);
    // Copy result
    memcpy(f, fval, *m * sizeof(double));
    // Clean up
    mxDestroyArray(mlData->plhs[0]);
}



//
// Check User Arguments
//
void MKLJac::checkInputArgs(const opti_mex_utils::OptiMexArgs& args, bool& haveSize, bool& haveTol)
{   
    // Defaults
    haveTol  = false;
    haveSize = false;

    if (args.nrhs < 2)
    {
        OPTIMex::error("You must supply at least 2 arguments to mklJac");        
    }
    if (OPTIMex::isFunction(pFUN) == false)
    {
        OPTIMex::error("The first argument (fun) must be a MATLAB function handle");
    }
    if ((OPTIMex::isDoubleVector(pX) == false) && (OPTIMex::isDoubleScalar(pX) == false))
    {
        OPTIMex::error("The second argument (x) must be a real double vector");
    }
    if (args.nrhs > 2)
    {
        if (OPTIMex::isDoubleScalar(pLEN) == false)
        {
            OPTIMex::error("The third argument (numRow), if specified, must be a double scalar");
        }
        double x = OPTIMex::getDoubleScalar(pLEN);
        if ((x < 1.0) || (x > 1e8))
        {
            OPTIMex::error("The third argument (numRow), if specified, must be 1 <= numRow <= 1e8");
        }
        haveSize = true;
        if (args.nrhs > 3)
        {
            if (OPTIMex::isDoubleScalar(pTOL) == false)
            {
                OPTIMex::error("The fourth argument (tol), if specified, must be a double scalar");
            }
            double tol = OPTIMex::getDoubleScalar(pTOL);
            if ((tol < std::numeric_limits<double>::epsilon()) || (tol > 1.0))
            {
                OPTIMex::error("The fourth argument (tol), if specified, must be eps(1) <= tol <= 1");
            }
            haveTol = true;
        }
    }
}