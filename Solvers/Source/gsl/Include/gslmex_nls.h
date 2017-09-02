/* GSLMEX_NLS - A MATLAB MEX Interface to GSL NLS Solver
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "mex.h"
#include "gslmex.h"
#include "opti_mex_utils.h"
#include <gsl_multifit_nlinear.h>

namespace opti_gsl_nls
{

// A Class for Solving NLS problems using GSL
class GSLMexNLS
{
    public:  
        static void solve(const opti_mex_utils::OptiMexArgs& args);

    private:
        static void checkInputArgs(const opti_mex_utils::OptiMexArgs& args);
        static void printIter(const gsl_multifit_nlinear_workspace* gslWs, double execTime);
        static void printHeader(const opti_mex_utils::OptiSolverProperties& solverInfo, const gsl_multifit_nlinear_workspace* gslWs);
        static void printFooter(int status, int info, double fval, double niter, double exitflag);
};

} // namespace opti_gsl_nls