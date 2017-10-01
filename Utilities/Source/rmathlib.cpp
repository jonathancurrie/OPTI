/* RMATHLIBMEX - A MATLAB MEX Interface to the R Math Library
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2012-2017
 * www.inverseproblem.co.nz
 */
#include "rmathlib.h"
#include <algorithm>
using namespace opti_mex_utils;
            

// Solver Defines
#define SOLVER_NAME     ("rmathlib: R Math Library")
#define SOLVER_LICENSE  (OptiSolverLicense::GPL)
#define SOLVER_LINK     ("http://cran.stat.sfu.ca/")
#define SOLVER_VERSION  (R_VERSION_STRING)

//
// Main Function
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    OptiSolverProperties solver(SOLVER_NAME, SOLVER_VERSION, SOLVER_LINK, SOLVER_LICENSE);
    OptiMexArgs args(solver, prhs, nrhs, plhs, nlhs);

    // Check if just requesting version/solver info
    if (OPTIMex::checkNoArgs(args) == true)
    {
        return; 
    }

    // Perform the calculation and return
    opti_utility::RMathLib rmath(args);
}  


namespace opti_utility
{
// Inputs/Outputs
#define pFUN    (args.prhs[0])
#define pIN1    (args.prhs[1])
#define pIN2    (args.prhs[2])
#define pIN3    (args.prhs[3])
#define pIN4    (args.prhs[4])
#define pIN5    (args.prhs[5])
#define pSOL    (args.plhs[0])


//
// Constructor (builds + solves)
//
RMathLib::RMathLib(const OptiMexArgs& args)
{    
    // Check input args
    RFcnType fcnType    = RFcnType::D1;
    std::string fcnName = checkInputArgs(args, fcnType);

    // Create Outputs
    pSOL = OPTIMex::createDoubleMatrix(mNumRows, mNumCols);
    double* v = mxGetPr(pSOL);
    
    // Solve based on supplied args
    switch (fcnType)
    {
        case RFcnType::D1:   { solve(mD1Fcns[fcnName], pIN1, v); break; }
        case RFcnType::D2:   { solve(mD2Fcns[fcnName], pIN1, pIN2, v); break; }
        case RFcnType::D2I1: { solve(mD2I1Fcns[fcnName], pIN1, pIN2, v); break; }
        case RFcnType::D2I2: { solve(mD2I2Fcns[fcnName], pIN1, pIN2, v); break; }
        case RFcnType::D3I1: { solve(mD3I1Fcns[fcnName], pIN1, pIN2, pIN3, v); break; }
        case RFcnType::D3I2: { solve(mD3I2Fcns[fcnName], pIN1, pIN2, pIN3, v); break; }
    }
}

//
// Specific Solving Routines
//
// D1
void RMathLib::solve(std::function<double(double)> fcn, const mxArray* in1, double* v) const
{
    double x = OPTIMex::getDoubleScalar(in1);

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            v[r + c*mNumRows] = fcn(x);
        }
    }
}

// D2
void RMathLib::solve(std::function<double(double,double)> fcn, const mxArray* in1, const mxArray* in2, double* v) const
{
    double x1 = OPTIMex::getDoubleScalar(in1);
    double x2 = OPTIMex::getDoubleScalar(in2);

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            v[r + c*mNumRows] = fcn(x1, x2);
        }
    }
}

// D2I1
void RMathLib::solve(std::function<double(double,double,int)> fcn, const mxArray* in1, const mxArray* in2, double* v) const
{
    double* x1 = mxGetPr(in1);
    double* x2 = mxGetPr(in2);
    size_t  x1mul = 1;
    size_t  x2mul = 1;

    // Determine index adjustment based on scalar expansion
    switch (mScalarExpReq)
    {
        case ScalarExp::NO_EXP:     { break; }
        case ScalarExp::IN1_SCALAR: { x1mul = 0; break; }
        case ScalarExp::IN2_SCALAR: { x2mul = 0; break; }
        default: { OPTIMex::error("Scalar expansion method not supported"); break; }
    }

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            size_t idx = r + c*mNumRows;
            v[idx] = fcn(x1[idx*x1mul], x2[idx*x2mul], mLogP);
        }
    }
}

// D2I2
void RMathLib::solve(std::function<double(double,double,int,int)> fcn, const mxArray* in1, const mxArray* in2, double* v) const
{
    double* x1 = mxGetPr(in1);
    double* x2 = mxGetPr(in2);
    size_t  x1mul = 1;
    size_t  x2mul = 1;

    // Determine index adjustment based on scalar expansion
    switch (mScalarExpReq)
    {
        case ScalarExp::NO_EXP:     { break; }
        case ScalarExp::IN1_SCALAR: { x1mul = 0; break; }
        case ScalarExp::IN2_SCALAR: { x2mul = 0; break; }
        default: { OPTIMex::error("Scalar expansion method not supported"); break; }
    }

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            size_t idx = r + c*mNumRows;
            v[idx] = fcn(x1[idx*x1mul], x2[idx*x2mul], mLowerTail, mLogP);
        }
    }
}

// D3I1
void RMathLib::solve(std::function<double(double,double,double,int)> fcn, const mxArray* in1, const mxArray* in2, const mxArray* in3, double* v) const
{
    double* x1 = mxGetPr(in1);
    double* x2 = mxGetPr(in2);
    double* x3 = mxGetPr(in3);
    size_t  x1mul = 1;
    size_t  x2mul = 1;
    size_t  x3mul = 1;

    // Determine index adjustment based on scalar expansion
    switch (mScalarExpReq)
    {
        case ScalarExp::NO_EXP:     { break; }
        case ScalarExp::IN1_SCALAR: { x1mul = 0; break; }
        case ScalarExp::IN2_SCALAR: { x2mul = 0; break; }
        case ScalarExp::IN3_SCALAR: { x3mul = 0; break; }
        case ScalarExp::IN12_SCALAR: { x1mul = 0; x2mul = 0; break; }
        case ScalarExp::IN13_SCALAR: { x1mul = 0; x3mul = 0; break; }
        case ScalarExp::IN23_SCALAR: { x2mul = 0; x3mul = 0; break; }
        default: { OPTIMex::error("Scalar expansion method not supported"); break; }
    }

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            size_t idx = r + c*mNumRows;
            v[idx] = fcn(x1[idx*x1mul], x2[idx*x2mul], x3[idx*x3mul], mLogP);
        }
    }
}

// D3I2
void RMathLib::solve(std::function<double(double,double,double,int,int)> fcn, const mxArray* in1, const mxArray* in2, const mxArray* in3, double* v) const
{
    double* x1 = mxGetPr(in1);
    double* x2 = mxGetPr(in2);
    double* x3 = mxGetPr(in3);
    size_t  x1mul = 1;
    size_t  x2mul = 1;
    size_t  x3mul = 1;

    // Determine index adjustment based on scalar expansion
    switch (mScalarExpReq)
    {
        case ScalarExp::NO_EXP:     { break; }
        case ScalarExp::IN1_SCALAR: { x1mul = 0; break; }
        case ScalarExp::IN2_SCALAR: { x2mul = 0; break; }
        case ScalarExp::IN3_SCALAR: { x3mul = 0; break; }
        case ScalarExp::IN12_SCALAR: { x1mul = 0; x2mul = 0; break; }
        case ScalarExp::IN13_SCALAR: { x1mul = 0; x3mul = 0; break; }
        case ScalarExp::IN23_SCALAR: { x2mul = 0; x3mul = 0; break; }
        default: { OPTIMex::error("Scalar expansion method not supported"); break; }
    }

    // Loop and calculate
    for (size_t r = 0; r < mNumRows; r++)
    {
        for (size_t c = 0; c < mNumCols; c++)
        {
            size_t idx = r + c*mNumRows;
            v[idx] = fcn(x1[idx*x1mul], x2[idx*x2mul], x3[idx*x3mul], mLowerTail, mLogP);
        }
    }
}


//
// Check Input Arguments
//
std::string RMathLib::checkInputArgs(const opti_mex_utils::OptiMexArgs& args, RFcnType& fcnType)
{    
    // Check common inputs
    OPTIMex::checkNumArgsIn(args.nrhs, 3, "rmathlib", "rmathlib('fcn',arg1,arg2)");
    OPTIMex::checkIsString(pFUN, "Function");
    OPTIMex::checkIsDouble(pIN1, "Argument 1");
    OPTIMex::checkIsDouble(pIN2, "Argument 2");
    
    // Collect the Function Name
    std::string fcnName = OPTIMex::getString(pFUN);
    // Match to Function Type
    fcnType = matchFcn(fcnName);
    
    // Depending on function found, check for correct number of args and sizes
    switch (fcnType)
    {
        case RFcnType::D1:
        {
            OPTIMex::checkNumArgsIn(args.nrhs, 4, "rmathlib", "rmathlib('fun',df,m,n)");
            OPTIMex::checkIsDoubleScalar(pIN1, "Argument 1 (degrees of freedom)");
            OPTIMex::checkIsDoubleScalarInBounds(pIN2, 1.0, 1e12, "Argument 2 (# rows)");
            OPTIMex::checkIsDoubleScalarInBounds(pIN3, 1.0, 1e12, "Argument 3 (# cols)");
            // Update Sizes
            mNumRows = static_cast<size_t>(OPTIMex::getDoubleScalar(pIN2));
            mNumCols = static_cast<size_t>(OPTIMex::getDoubleScalar(pIN3));            
            break;
        }
        case RFcnType::D2:
        {
            OPTIMex::checkNumArgsIn(args.nrhs, 5, "rmathlib", "rmathlib('fun',df1,df2,m,n)");
            OPTIMex::checkIsDoubleScalar(pIN1, "Argument 1 (degrees of freedom 1)");
            OPTIMex::checkIsDoubleScalar(pIN2, "Argument 2 (degrees of freedom 2)");
            OPTIMex::checkIsDoubleScalarInBounds(pIN3, 1.0, 1e12, "Argument 3 (# rows)");
            OPTIMex::checkIsDoubleScalarInBounds(pIN4, 1.0, 1e12, "Argument 4 (# cols)");
            // Update Sizes
            mNumRows = static_cast<size_t>(OPTIMex::getDoubleScalar(pIN3));
            mNumCols = static_cast<size_t>(OPTIMex::getDoubleScalar(pIN4));            
            break;
        }
        case RFcnType::D2I1:
        {
            OPTIMex::checkIsDouble(pIN1, "Argument 1");
            OPTIMex::checkIsDouble(pIN2, "Argument 2");
            if ((args.nrhs > 3) && (OPTIMex::isScalar(pIN3) == true))
            {
                mLogP = OPTIMex::getLogicalScalar(pIN3);
            }
            mScalarExpReq = OPTIMex::checkScalarExpansion(pIN1, pIN2, mNumRows, mNumCols);
            break;
        }
        case RFcnType::D2I2:
        {
            OPTIMex::checkIsDouble(pIN1, "Argument 1");
            OPTIMex::checkIsDouble(pIN2, "Argument 2");
            if ((args.nrhs > 3) && (OPTIMex::isScalar(pIN3) == true))
            {
                mLowerTail = OPTIMex::getLogicalScalar(pIN3);
            }
            if ((args.nrhs > 4) && (OPTIMex::isScalar(pIN4) == true))
            {
                mLogP = OPTIMex::getLogicalScalar(pIN4);
            }
            mScalarExpReq = OPTIMex::checkScalarExpansion(pIN1, pIN2, mNumRows, mNumCols);
            break;            
        }
        case RFcnType::D3I1:
        {
            OPTIMex::checkNumArgsIn(args.nrhs, 4, "rmathlib", "rmathlib('fun',x,df1,df2)");
            OPTIMex::checkIsDouble(pIN1, "Argument 1");
            OPTIMex::checkIsDouble(pIN2, "Argument 2");
            OPTIMex::checkIsDouble(pIN3, "Argument 3");
            if ((args.nrhs > 4) && (OPTIMex::isScalar(pIN4) == true))
            {
                mLogP = OPTIMex::getLogicalScalar(pIN4);
            }
            mScalarExpReq = OPTIMex::checkScalarExpansion(pIN1, pIN2, pIN3, mNumRows, mNumCols);
            break; 
        }
        case RFcnType::D3I2:
        {
            OPTIMex::checkNumArgsIn(args.nrhs, 4, "rmathlib", "rmathlib('fun',x,df1,df2)");
            OPTIMex::checkIsDouble(pIN1, "Argument 1");
            OPTIMex::checkIsDouble(pIN2, "Argument 2");
            OPTIMex::checkIsDouble(pIN3, "Argument 3");
            if ((args.nrhs > 4) && (OPTIMex::isScalar(pIN4) == true))
            {
                mLowerTail = OPTIMex::getLogicalScalar(pIN4);
            }
            if ((args.nrhs > 5) && (OPTIMex::isScalar(pIN5) == true))
            {
                mLogP = OPTIMex::getLogicalScalar(pIN5);
            }
            mScalarExpReq = OPTIMex::checkScalarExpansion(pIN1, pIN2, pIN3, mNumRows, mNumCols);
            break; 
        }
    }
    // Return the loaded function name
    return fcnName;
}

//
// Match Name to RMathLib Function Type
//
RFcnType RMathLib::matchFcn(std::string& fcnName) const
{
    // Convert the string to lowercase
    std::transform(fcnName.begin(), fcnName.end(), fcnName.begin(), ::tolower);

    if (mD1Fcns.count(fcnName) > 0)
    {
        return RFcnType::D1;
    }
    else if (mD2Fcns.count(fcnName) > 0)
    {
        return RFcnType::D2;
    }
    else if (mD2I1Fcns.count(fcnName) > 0)
    {
        return RFcnType::D2I1;
    }
    else if (mD2I2Fcns.count(fcnName) > 0)
    {
        return RFcnType::D2I2;
    }
    else if (mD3I1Fcns.count(fcnName) > 0)
    {
        return RFcnType::D3I1;
    }
    else if (mD3I2Fcns.count(fcnName) > 0)
    {
        return RFcnType::D3I2;
    }
    else
    {
        OPTIMex::error("OPTIMex:InputError", "Unknown RMathLib Function '%s'", fcnName.c_str());
        return RFcnType::D1; // not used 
    }
}

} // namespace opti_utility
