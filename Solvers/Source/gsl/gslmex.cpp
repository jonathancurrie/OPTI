/* GSLMEX - A MATLAB MEX Interface to GSL
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "gslmex.h"
#include "opti_util.h"
#ifdef LINK_MKL
#include <mkl.h>
#endif

// Local Function Prototypes
void printSolverInfo();
GslProbType readProblemType(const mxArray *prhs[], int nrhs);

// Main Function
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    // No input macro
    if(nrhs < 1) {
        if(nlhs < 1)
            printSolverInfo();
        else
        {
            plhs[0] = mxCreateString(GSL_VERSION);
            plhs[1] = mxCreateDoubleScalar(OPTI_VER);
        }
        return;
    }    
    
    // Get the problem type we are solving
    GslProbType prob = readProblemType(prhs, nrhs);
    
    // Solve based on the problem type
    switch (prob)
    {
        // Nonlinear least squares
        case GslProbType::NLS:
        {
            gslSolveNLS(prhs, nrhs, plhs, nlhs);
            break;
        }
    }        
}             

// Determine Problem Type 
GslProbType readProblemType(const mxArray *prhs[], int nrhs)
{
    if (mxIsStruct(prhs[0]) == false)
    {
        mexErrMsgTxt("The first argument to gsl() must be a structure");
        return GslProbType::UNKNOWN;
    }
    mxArray* probType = mxGetField(prhs[0], 0, "probType");
    
    if ((probType != nullptr) && !mxIsEmpty(probType) && mxIsChar(probType))
    {
        char probTypeStr[256];
        mxGetString(probType, probTypeStr, 256);
        if (_stricmp(probTypeStr,"NLS") == 0)
        {
            return GslProbType::NLS;
        }
        else
        {
            MEX_ERR("Unknown problem type '%s'", probTypeStr);
        }
    }
    else
    {
        mexErrMsgTxt("The problem structure must contain the field 'probType'");
    }
    
    return GslProbType::UNKNOWN;
}

//Print Solver Information
void printSolverInfo()
{     
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" GSL: GNU Scientific Library [v%s]\n", GSL_VERSION);
    PRINT_BUILD_INFO;
    mexPrintf("  - Released under the GNU General Public License: https://www.gnu.org/copyleft/gpl.html\n");
    mexPrintf("  - Source available from: https://www.gnu.org/software/gsl/\n\n");
    
    mexPrintf(" This binary is statically linked to the following software:\n");
    #ifdef LINK_NETLIB_BLAS
        mexPrintf("  - NETLIB BLAS: http://www.netlib.org/blas/\n  - NETLIB LAPACK: http://www.netlib.org/lapack/\n");
    #endif
    #ifdef LINK_MKL
        mexPrintf("  - Intel Math Kernel Library [v%d.%d R%d]\n",__INTEL_MKL__,__INTEL_MKL_MINOR__,__INTEL_MKL_UPDATE__);        
    #endif

    mexPrintf("\n MEX Interface J.Currie 2017 [BSD3] (www.inverseproblem.co.nz)\n");
    mexPrintf("-----------------------------------------------------------\n");
}