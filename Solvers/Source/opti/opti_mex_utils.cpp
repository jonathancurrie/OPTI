/* OPTI MEX Utilities Source File
 * (C) Inverse Problem Limited 2017
 * J. Currie 
 */

#include "opti_mex_utils.h"
#include <stdarg.h>

//
// Extern Defines for Ctrl-C Detection (only supported for MSVC compiler)
//
#ifdef _MSC_VER
// Ctrl-C Detection 
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
#endif

// External Libs Headers as Included
#ifdef LINK_MKL
#include "mkl.h"
#endif

namespace opti_mex_utils
{
//
// Data Creation
//
// Create Double Precision Scalar
mxArray* OPTIMex::createDoubleScalar(double val)
{
    return mxCreateDoubleScalar(val);
}

// Create Double Precision Matrix
mxArray* OPTIMex::createDoubleMatrix(size_t nrow, size_t ncol)
{
    return mxCreateDoubleMatrix(nrow, ncol, mxREAL);
}

// Create Structure with Specified Fields
mxArray* OPTIMex::createStruct(const char* fieldNames[], size_t numFields)
{
    return mxCreateStructMatrix(1, 1, numFields, fieldNames);
}

// Create Field in Structure - Scalar
double* OPTIMex::createFieldScalar(mxArray* data, const char* fieldName, double val)
{
    if (OPTIMex::isValidStruct(data))
    { 
        mxArray* field = OPTIMex::createDoubleScalar(val);
        mxSetField(data, 0, fieldName, field);
        return mxGetPr(field);
    }
    else
    {
        OPTIMex::error("Cannot create double scalar within structure field '%s'", fieldName);
    }
}

// Create Field in Structure - Matrix
double* OPTIMex::createFieldMatrix(mxArray* data, const char* fieldName, size_t nrow, size_t ncol)
{
    if (OPTIMex::isValidStruct(data))
    { 
        mxArray* field = OPTIMex::createDoubleMatrix(nrow, ncol);
        mxSetField(data, 0, fieldName, field);
        return mxGetPr(field);
    }
    else
    {
        OPTIMex::error("Cannot create double matrix within structure field '%s'", fieldName);
    }
}


//
// Data Access
//
// Return how many elements in specified structure field
size_t OPTIMex::getFieldNumElem(const mxArray* data, const char* field)
{
    return getNumElem(getField(data, field));
}

size_t OPTIMex::getNumElem(const mxArray* data)
{
    if (isEmpty(data) == false)
    {
        return mxGetNumberOfElements(data); 
    }
    return 0;
}

size_t OPTIMex::getNumRows(const mxArray* data)
{
    if (isEmpty(data) == false)
    {
        return mxGetN(data); 
    }
    return 0;
}

size_t OPTIMex::getNumCols(const mxArray* data)
{
    if (isEmpty(data) == false)
    {
        return mxGetM(data); 
    }
    return 0;
}

// Access Structure Field
mxArray* OPTIMex::getField(const mxArray* data, const char* fieldName)
{
    if (OPTIMex::isValidField(data, fieldName))
    {
        return mxGetField(data, 0, fieldName);
    }
    else
    {
        OPTIMex::error("Cannot access field '%s' in structure", fieldName);
    }
}

// Access double array in Structure Field
double* OPTIMex::getFieldPr(const mxArray* data, const char* fieldName)
{
    mxArray* field = getField(data, fieldName);
    if (isRealDouble(field))
    {
        return mxGetPr(field);
    }
    else
    {
        OPTIMex::error("Field '%s' does not contain a double variable", fieldName);
    }
}



//
// Options Access
//
int OPTIMex::getIntegerOption(const mxArray* opts, const char* optionName, int& option)
{
    if (OPTIMex::isValidField(opts, optionName))
    { 
        option = (int)*mxGetPr(mxGetField(opts, 0, optionName));        
        return OPTI_MEX_SUCCESS;
    }
    return OPTI_MEX_FAILURE;
}
int OPTIMex::getIntegerOption(const mxArray* opts, const char* optionName, size_t& option)
{
    if (OPTIMex::isValidField(opts, optionName))
    { 
        option = (size_t)*mxGetPr(mxGetField(opts, 0, optionName));        
        return OPTI_MEX_SUCCESS;
    }
    return OPTI_MEX_FAILURE;
}
int OPTIMex::getDoubleOption(const mxArray* opts, const char* optionName, double& option)
{
    if (OPTIMex::isValidField(opts, optionName))
    { 
        option = *mxGetPr(mxGetField(opts, 0, optionName));        
        return OPTI_MEX_SUCCESS;
    }
    return OPTI_MEX_FAILURE;
}


//
// Data Validation
//
bool OPTIMex::isValidStruct(const mxArray* data)
{
    if (data != nullptr)
    {
        return mxIsStruct(data) && !mxIsEmpty(data);
    }
    else
    {
        return false;
    }
}

bool OPTIMex::isValidField(const mxArray* data, const char* field)
{
    if (OPTIMex::isValidStruct(data) && (field != nullptr))
    {
        return mxGetField(data, 0, field) != nullptr;
    }
    else
    {
        return false;
    }
}

bool OPTIMex::isString(const mxArray* data)
{
    if (data != nullptr)
    {
        return mxIsChar(data) && !mxIsEmpty(data);
    }
    else
    {
        return false;
    }
}

bool OPTIMex::isDoubleScalar(const mxArray* data)
{
    if (isRealDouble(data) && (getNumElem(data) == 1))
    {
        return true;
    }
    return false;
}

bool OPTIMex::isDoubleVector(const mxArray* data)
{
    if (isRealDouble(data) && isVector(data))
    {
        return true;
    }
    return false;
}

bool OPTIMex::isDoubleMatrix(const mxArray* data)
{
    if (isRealDouble(data) && isMatrix(data))
    {
        return true;
    }
    return false;
}

bool OPTIMex::isRealDouble(const mxArray* data)
{
    if (isEmpty(data) == false)
    {
        return (mxIsComplex(data) == false) && (mxGetClassID(data) == mxDOUBLE_CLASS);
    }
    return false;
}

bool OPTIMex::isEmpty(const mxArray* data)
{
    return mxIsEmpty(data);
}

bool OPTIMex::isMatrix(const mxArray* data)
{
    return (getNumRows(data) > 0) && (getNumCols(data) > 0);
}

bool OPTIMex::isVector(const mxArray* data)
{
    return (getNumRows(data) > 0) ^ (getNumCols(data) > 0);
}

bool OPTIMex::isScalar(const mxArray* data)
{
    return getNumElem(data) == 1;
}



//
// Error Message
//
void OPTIMex::error(const char* format, ...)
{
    char errBuf[1024];
    va_list args;
    va_start(args, format);
    vsnprintf(errBuf, 1024, format, args);
    va_end(args);
    
    // Generates Exception
    mexErrMsgTxt(errBuf);
}

//
// Check if Ctrl-C Has Been Pressed
//
bool OPTIMex::ctrlCPressed(const char* solverName)
{
    if (utIsInterruptPending())
    {
        utSetInterruptPending(false); // Clear Ctrl-C Status
        mexPrintf("\nCtrl-C Detected. Exiting %s...\n\n", solverName);
        return true;
    }
    return false;
}


//
// Check OPTI Termination Checks
//
int OPTIMex::checkTermination(OPTITimer& timer, size_t numIter, size_t numFeval, size_t numNodes, const OPTITermSettings& termSettings)
{
    if (timer.execTime() > termSettings.maxTime)
    {
        return OPTI_TERM_MAX_TIME;
    }

    if (numFeval > termSettings.maxFeval)
    {
        return OPTI_TERM_MAX_FEVAL;
    }

    if (numIter > termSettings.maxIter)
    {
        return OPTI_TERM_MAX_ITER;
    }

    if (numNodes > termSettings.maxNodes)
    {
        return OPTI_TERM_MAX_NODES;
    }

    if (ctrlCPressed(termSettings.solverName.c_str()) == true)
    {
        return OPTI_TERM_USER_EXIT;
    }

    return OPTI_TERM_CONTINUE;
}

//
// Callback Memory Initialization
//
void OPTIMex::initMatlabCallbackData(matlab_cb_data_t& callbackData, const mxArray* problemData, size_t ndec, size_t ndata)
{
    // Get the function - will throw an exception if it doesn't exist
    callbackData.prhs_f[0] = getField(problemData, "fun");
    callbackData.prhs_f[1] = createDoubleMatrix(ndec, 1); // x0

    // If we have a gradient, get it
    if (isValidField(problemData, "grad"))
    {
        callbackData.prhs_g[0] = getField(problemData, "grad");
        callbackData.haveGrad  = true;
        callbackData.prhs_g[1] = createDoubleMatrix(ndec, 1); // x0
    }

    // Assign data fitting values, if present
    if (isValidField(problemData, "ydata"))
    {
        callbackData.ydata = mxGetPr(getField(problemData, "ydata"));
    }

    // Assign Sizes
    callbackData.ndec    = ndec;
    callbackData.ndata   = ndata;
}

void OPTIMex::initIterFunCallbackData(iterfun_cb_data_t& callbackData, mxArray* iterFunHandle, size_t ndec)
{
    callbackData.prhs[0] = iterFunHandle;
    callbackData.prhs[1] = createDoubleScalar(0.0);     // niter
    callbackData.prhs[2] = createDoubleScalar(0.0);     // fval
    callbackData.prhs[3] = createDoubleMatrix(ndec, 1); // x0
}


//
// MATLAB Callbacks
//
double* OPTIMex::callMatlabObjective(matlab_cb_data_t* mlData, double* x)
{
    // Copy in x
    mlData->plhs[0] = nullptr;
    memcpy(mxGetPr(mlData->prhs_f[1]), x, mlData->ndec * sizeof(double));

    // Call MATLAB
    if(mexCallMATLAB(1, mlData->plhs, 2, mlData->prhs_f, mlData->f) != 0)
    {
        mexErrMsgTxt("Error calling Objective Function!");
    }
    mlData->nfeval++;
    
    // Access Objective Pointer
    return mxGetPr(mlData->plhs[0]);
}

double* OPTIMex::callMatlabGradient(matlab_cb_data_t* mlData, double* x)
{
    // Copy in x
    mlData->plhs[0] = nullptr;
    memcpy(mxGetPr(mlData->prhs_g[1]), x, mlData->ndec * sizeof(double));

    // Call MATLAB
    if(mexCallMATLAB(1, mlData->plhs, 2, mlData->prhs_g, mlData->f) != 0)
    {
        mexErrMsgTxt("Error calling Gradient Function!");
    }
    mlData->ngeval++;
    
    // Access Gradient Pointer
    return mxGetPr(mlData->plhs[0]);
}

bool OPTIMex::callMatlabIterFun(iterfun_cb_data_t* iterData, size_t niter, double fval, double* x)
{        
    // Assign iteration data
    iterData->plhs[0] = NULL;        
    memcpy(mxGetData(iterData->prhs[1]), &niter, sizeof(int));
    memcpy(mxGetPr(iterData->prhs[2]), &fval, sizeof(double));
    memcpy(mxGetPr(iterData->prhs[3]), x, iterData->ndec * sizeof(double));
    
    // Call MATLAB
    if (mexCallMATLAB(1, iterData->plhs, 4, iterData->prhs, iterData->f) != 0)
    {
        mexErrMsgTxt("Error calling Iteration Callback Function!");
    }

    //Collect return argument
    bool stop = *(bool*)mxGetData(iterData->plhs[0]);

    // Clean up memory & return stop status
    mxDestroyArray(iterData->plhs[0]);
    return stop;
}


//
// Copy Helpers
//
void OPTIMex::copyTranspose(double* dest, double* src, size_t nrows, size_t ncols)
{
    for(size_t r = 0; r < nrows; ++r)
    {
        for(size_t c = 0; c < ncols; ++c)
        {
            dest[c + r*ncols] = src[r + c*nrows];
        }        
    }    
}


//
// Check for no user arguments, return version info / print solver info as required
//
bool OPTIMex::checkNoArgs(const OptiMexArgs& args)
{
    if(args.nrhs < 1) 
    {
        if(args.nlhs < 1)
        {
            printSolverInfo(args.solverInfo);
        }
        else
        {
            args.plhs[0] = mxCreateString(args.solverInfo.solverVersion.c_str());
            #ifdef OPTI_VER
                args.plhs[1] = mxCreateDoubleScalar(OPTI_VER);
            #else
                args.plhs[1] = mxCreateDoubleScalar(NAN);
            #endif
        }
        return true;
    }  
    return false;
}

//
// Print Solver Header
//
void OPTIMex::printSolverHeader(const OptiSolverProperties& solverInfo, const char* algorithmString, size_t ndec, size_t ndata)
{
    mexPrintf("\n------------------------------------------------------------------\n");        
    mexPrintf(" This is %s v%s\n", solverInfo.solverName.c_str(), solverInfo.solverVersion.c_str());           
    mexPrintf(" Algorithm: %s", algorithmString);
    mexPrintf(" Problem Properties:\n");
    mexPrintf(" # Decision Variables:     %4d\n", ndec);
    mexPrintf(" # Data Points:            %4d\n", ndata);
    mexPrintf("------------------------------------------------------------------\n");
}

//
// Solver Iteration Print Out
//
void OPTIMex::printSolverIter(size_t niter, double execTime, const char *header, const char *format, ...)
{
    char msgBuf[1024];
    // Display header every n iters
    if(!(niter % OPTI_MEX_HEADER_PRINT_FREQ) || (niter == 1))
    {
        // Generate Full Header
        strcpy(msgBuf, " iter              time            ");
        strcat(msgBuf, header);
        mexPrintf(msgBuf);
        mexEvalString("drawnow;"); // ensure flushed every n iters
    }        
    
    // Generate Start of detail line
    int n = snprintf(msgBuf, 1024, "%5d         %9.3g   " , niter, execTime);

    // Add rest of the detail
    va_list args;
    va_start(args, format);
    vsnprintf(&msgBuf[n], 1024-n, format, args);
    va_end(args);

    // Print it 
    mexPrintf(msgBuf);    
}

//
// Print Solver Information
//
void OPTIMex::printSolverInfo(const OptiSolverProperties& solver)
{
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" %s [v%s]\n", solver.solverName.c_str(), solver.solverVersion.c_str());
    PRINT_BUILD_INFO;
    std::string licenseStr;
    switch (solver.solverLicense)
    {
        case OptiSolverLicense::GPL: { licenseStr = OPTI_GPL_LICENSE_LINK; break; }
        case OptiSolverLicense::EPL: { licenseStr = OPTI_EPL_LICENSE_LINK; break; }
        case OptiSolverLicense::LGPL: { licenseStr = OPTI_LGPL_LICENSE_LINK; break; }
    }
    if (licenseStr.empty() == false)
    {
        mexPrintf("  - Released under the %s\n", licenseStr.c_str());
    }
    mexPrintf("  - Source available from: %s\n\n", solver.solverLink.c_str());
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef SOPLEX_VERSION
        mexPrintf("  - SoPlex [v%d] (ZIB Academic License)\n",SOPLEX_VERSION);
    #endif
    #ifdef IPOPT_VERSION
        mexPrintf("  - Ipopt  [v%s] (Eclipse Public License)\n",IPOPT_VERSION);
    #endif
    #ifdef CPPAD_PACKAGE_STRING
        strcpy(msgbuf,&CPPAD_PACKAGE_STRING[6]);
        mexPrintf("  - CppAD  [v%s] (Eclipse Public License)\n",msgbuf); 
    #endif
    #ifdef LINK_MUMPS
        mexPrintf("  - MUMPS  [v%s]\n",MUMPS_VERSION);
        mexPrintf("  - METIS  [v4.0.3] (Copyright University of Minnesota)\n");
    #endif
    #ifdef LINK_ASL
        mexPrintf("  - ASL    [v%d] (Netlib)\n",ASLdate_ASL);
    #endif
    #ifdef LINK_NETLIB_BLAS
        mexPrintf("  - NETLIB BLAS: http://www.netlib.org/blas/\n  - NETLIB LAPACK: http://www.netlib.org/lapack/\n");
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
        #if defined(LINK_METIS) && !defined(LINK_MUMPS)
            mexPrintf("  - MeTiS [v4.0.3] Copyright University of Minnesota\n");
        #endif
    #endif    
    #if defined(LINK_ML_MA57) || defined(LINK_PARDISO) || defined(LINK_CPLEX)
        mexPrintf("\n And is dynamically linked to the following software:\n");
        #ifdef LINK_ML_MA57
            mexPrintf("  - MA57    [v3.0] (Included as part of the MATLAB distribution)\n");
        #endif
        #ifdef LINK_PARDISO //BASEL VERSION
            mexPrintf("  - PARDISO [v4.1.2]\n");
        #endif
        #ifdef LINK_CPLEX
            mexPrintf("  - CPLEX  [v%d.%d.%d]\n",CPX_VERSION_VERSION,CPX_VERSION_RELEASE,CPX_VERSION_MODIFICATION);
        #endif
    #endif
            
    #ifdef HAVE_LINEARSOLVERLOADER
        mexPrintf("\n Dynamically Linked to libhsl.dll containing the HSL Numerical Routines (MA27, MA77, etc)\n   (NOTE: This must exist on the MATLAB path)\n");
    #endif

    mexPrintf("\n MEX Interface J.Currie 2012-2017 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}


//
// OPTI Timer
//
OPTITimer::OPTITimer(bool startNow)
{
    if (startNow == true)
    {
        start();
    }
    else
    {
        _counterPeriod       = 0.0;
        _startCount.QuadPart = 0;
        _configured          = false;
    }
}

bool OPTITimer::start(void)
{
    LARGE_INTEGER freq;
    if ((_configured = QueryPerformanceFrequency(&freq) && QueryPerformanceCounter(&_startCount)))
    {
        _counterPeriod = 1.0 / static_cast<double>(freq.QuadPart);
    }
    return _configured;
}

double OPTITimer::execTime(void)
{
    LARGE_INTEGER currentCount;
    if (_configured && QueryPerformanceCounter(&currentCount))
    {
        return static_cast<double>(currentCount.QuadPart - _startCount.QuadPart) * _counterPeriod;
    }
    else
    {
        return NAN;
    }
}

} // namespace opti_mex_utils