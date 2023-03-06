/* MKLJAC - A MATLAB MEX Interface to Intel MKL's DJACOBI
 * Released Under the BSD 3-Clause License:
 * https://www.controlengineering.co.nz/Wikis/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.controlengineering.co.nz
 */

#include "opti_mex_utils.h"

namespace opti_utility
{

class MKLJac
{
    public:  
        static void differentiate(const opti_mex_utils::OptiMexArgs& args);
    private:
        static void checkInputArgs(const opti_mex_utils::OptiMexArgs& args, bool& haveSize, bool& haveTol);
};

} // namespace opti_utility