/* GSLMEX - A MATLAB MEX Interface to GSL
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2017
 * www.inverseproblem.co.nz
 */

#pragma once

#include "opti_mex_utils.h"

namespace opti_gsl
{
// Common Return Codes
#define GSL_USER_EXIT   (-5)
#define GSL_MAX_TIME    (-6)
#define GSL_MAX_FEVAL   (-7)
#define GSL_NO_SOL      (-8)

// Problem Type Enum
enum class GslProbType {UNKNOWN,NLS};

class GSLMex
{
    public:
        GSLMex(const opti_mex_utils::OptiMexArgs& args);

    private:
        GslProbType _readProblemType(const opti_mex_utils::OptiMexArgs& args); 
};

} // namespace opti_gsl