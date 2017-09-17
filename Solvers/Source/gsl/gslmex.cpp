/* GSLMEX - A MATLAB MEX Interface to GSL
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "gslmex.h"
#include "gslmex_nls.h"
#include "opti_mex_utils.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_version.h>

using namespace opti_mex_utils;
using namespace opti_gsl;

// Solver Defines
#define SOLVER_NAME     ("GSL: GNU Scientific Library")
#define SOLVER_VERSION  (GSL_VERSION)
#define SOLVER_LICENSE  (OptiSolverLicense::GPL)
#define SOLVER_LINK     ("https://www.gnu.org/software/gsl/")

// Local Function Prototypes
void error_handler(const char * reason, const char * file, int line, int gsl_errno);

// Main Function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    OptiSolverProperties solver(SOLVER_NAME, SOLVER_VERSION, SOLVER_LINK, SOLVER_LICENSE);
    OptiMexArgs args(solver, prhs, nrhs, plhs, nlhs);

    // Check if just requesting version/solver info
    if (OPTIMex::checkNoArgs(args) == true)
    {
        return; 
    }

    // Build a GSLMex Object and solve the problem
    GSLMex gsl(args);
}          


//
// GSL Mex Constructor
//
GSLMex::GSLMex(const OptiMexArgs& args)
{
    // Get the problem type we are solving
    GslProbType prob = _readProblemType(args);
    
    // Set the error handler
    gsl_set_error_handler(error_handler);

    // Solve based on the problem type
    switch (prob)
    {
        // Nonlinear least squares
        case GslProbType::NLS:
        {
            opti_gsl_nls::GSLMexNLS::solve(args);
            break;
        }

        default:
        {
            OPTIMex::error("Unknown problem type to solve - please check your input arguments.");
        }
    } 
}


//
// Determine Problem Type 
//
GslProbType GSLMex::_readProblemType(const OptiMexArgs& args)
{
    if (OPTIMex::isValidField(args.prhs[0], "probType"))
    {
        mxArray* probType = OPTIMex::getField(args.prhs[0], "probType");
        if (OPTIMex::isString(probType))
        {
            char probTypeStr[256];
            mxGetString(probType, probTypeStr, 256);
            if (_stricmp(probTypeStr,"NLS") == 0)
            {
                return GslProbType::NLS;
            }
            else
            {
                OPTIMex::error("Unknown problem type '%s'", probTypeStr);
            }
        }
        else
        {
            OPTIMex::error("The probType field must be a string");
        }
    }
    else
    {
        OPTIMex::error("The problem structure must contain the field 'probType'");
    }
    return GslProbType::UNKNOWN;
}


// Local Error Handler
void error_handler(const char * reason, const char * file, int line, int gsl_errno)
{
    OPTIMex::error("Error in GSL \"%s\", line %d: %s (code %d)", file, line, reason, gsl_errno);
}
    
    
    