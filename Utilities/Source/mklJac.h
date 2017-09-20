/* MKLJAC - A MATLAB MEX Interface to Intel MKL's DJACOBI
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "opti_mex_utils.h"

class MKLJac
{
    public:  
        static void differentiate(const opti_mex_utils::OptiMexArgs& args);
    private:
        static void checkInputArgs(const opti_mex_utils::OptiMexArgs& args, bool& haveSize, bool& haveTol);
};