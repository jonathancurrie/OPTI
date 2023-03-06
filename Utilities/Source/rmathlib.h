/* RMATHLIB - A MATLAB MEX Interface to the R Math Library
 * Released Under the BSD 3-Clause License:
 * https://www.controlengineering.co.nz/Wikis/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */

#include "opti_mex_utils.h"
#include <Rmath.h>
#include <functional>
#include <map>

namespace opti_utility
{

// R Function Type
enum class RFcnType{D1,D2,D2I1,D2I2,D3I1,D3I2};

//
// The main RMathLib Class
//
class RMathLib
{
    public:  
        RMathLib(const opti_mex_utils::OptiMexArgs& args);
       
    private:        
        std::string checkInputArgs(const opti_mex_utils::OptiMexArgs& args, RFcnType& fcnType);
        RFcnType RMathLib::matchFcn(std::string& fcnName) const;
        void solve(std::function<double(double)>, const mxArray* in1, double* v) const;
        void solve(std::function<double(double,double)>, const mxArray* in1, const mxArray* in2, double* v) const;
        void solve(std::function<double(double,double,int)>, const mxArray* in1, const mxArray* in2, double* v) const;
        void solve(std::function<double(double,double,int,int)> fcn, const mxArray* in1, const mxArray* in2, double* v) const;
        void solve(std::function<double(double,double,double,int)> fcn, const mxArray* in1, const mxArray* in2, const mxArray* in3, double* v) const;
        void solve(std::function<double(double,double,double,int,int)> fcn, const mxArray* in1, const mxArray* in2, const mxArray* in3, double* v) const;

        // Internal Properties
        size_t mNumRows = 0;
        size_t mNumCols = 0;
        bool mLowerTail = true;    // Default is lower tail
        bool mLogP      = false;
        opti_mex_utils::ScalarExp mScalarExpReq = opti_mex_utils::ScalarExp::NO_EXP;   

        // Function Maps
        std::map<std::string, std::function<double(double)>> mD1Fcns {
            {"rt", &rt},
            {"rchisq", &rchisq},
            {"rexp", &rexp},
            {"rpois", &rpois}
        };     
        std::map<std::string, std::function<double(double,double)>> mD2Fcns {
            {"rf", &rf},
            {"rnorm", &rnorm},
            {"runif", &runif},
            {"rbeta", &rbeta},
            {"rcauchy", &rcauchy},
            {"rgamma", &rgamma},
            {"rlnorm", &rlnorm},
            {"rweibull", &rweibull}
        };
        std::map<std::string, std::function<double(double,double,int)>> mD2I1Fcns {
            {"dt", &dt},
            {"dchisq", &dchisq},
            {"dexp", &dexp},
            {"dpois", &dpois}
        };
        std::map<std::string, std::function<double(double,double,int,int)>> mD2I2Fcns {
            {"pt", &pt},
            {"qt", &qt},
            {"pchisq", &pchisq},
            {"qchisq", &qchisq},
            {"pexp", &pexp},
            {"qexp", &qexp},
            {"ppois", &ppois},
            {"qpois", &qpois}
        };
        std::map<std::string, std::function<double(double,double,double,int)>> mD3I1Fcns {
            {"df", &df},
            {"dnorm", &dnorm},
            {"dunif", &dunif},
            {"dbeta", &dbeta},
            {"dcauchy", &dcauchy},
            {"dgamma", &dgamma},
            {"dlnorm", &dlnorm},
            {"dweibull", &dweibull}
        };
        std::map<std::string, std::function<double(double,double,double,int,int)>> mD3I2Fcns {
            {"pf", &pf},
            {"qf", &qf},
            {"pnorm", &pnorm},
            {"qnorm", &qnorm},
            {"punif", &punif},
            {"qunif", &qunif},
            {"pbeta", &pbeta},
            {"qbeta", &qbeta},
            {"pcauchy", &pcauchy},
            {"qcauchy", &qcauchy},
            {"pgamma", &pgamma},
            {"qgamma", &qgamma},
            {"plnorm", &plnorm},
            {"qlnorm", &qlnorm},
            {"pweibull", &pweibull},
            {"qweibull", &qweibull},
        };
};

} // namespace opti_utility 
