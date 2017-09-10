/* OPTI MEX Utilities Header File
 * (C) Inverse Problem Limited 2017
 * J. Currie 
 */

#ifndef OPTI_MEX_UTILS
#define OPTI_MEX_UTILS

#include "mex.h"
#include "opti_build_utils.h"
#include <limits>
#include <string>
#include <Windows.h>

namespace opti_mex_utils
{
//
// OPTI Defaults
//
#define OPTI_DEFAULT_TOL            (1e-7)
#define OPTI_MAX_ITER               (1000)
#define OPTI_MAX_TIME               (1000)
#define OPTI_MAX_NODES              (1000)
#define OPTI_MAX_FEVAL              (10000)
#define OPTI_MEX_HEADER_PRINT_FREQ  (10)

// Common Verbosity Levels
#define OPTI_VERBOSE_OFF            (0)
#define OPTI_VERBOSE_MIN            (1)
#define OPTI_VERBOSE_ITER           (2)

//
// Return Codes
//
#define OPTI_MEX_SUCCESS        (0)
#define OPTI_MEX_FAILURE        (-1)

// Termination Check Return Codes
#define OPTI_TERM_CONTINUE      (0)
#define OPTI_TERM_USER_EXIT     (-51)
#define OPTI_TERM_MAX_FEVAL     (-52)
#define OPTI_TERM_MAX_ITER      (-53)
#define OPTI_TERM_MAX_TIME      (-54)
#define OPTI_TERM_MAX_NODES     (-55)

//
// MATLAB Callback Structures
//
#define OPTI_MEX_MAX_FCN_NAME  (128) 

// General MATLAB Callback Data
typedef struct matlab_cb_data_s
{
    char *f = "feval";
    mxArray *plhs[1];       // Return data
    mxArray *prhs_f[2];     // Function (fcnName, x)
    mxArray *prhs_g[2];     // Gradient (fcnName, x)
    size_t ndec         = 0;
    // Statistics
    bool haveGrad       = false;
    size_t nfeval       = 0;
    size_t ngeval       = 0;
    // Data fitting problems
    size_t ndata        = 0;
    double *ydata       = nullptr;
} matlab_cb_data_t;

// Iteration Function Callback Data
typedef struct iterfun_cb_data_s
{
    char *f = "feval";
    mxArray *plhs[1];   // Return data (stop logical)
    mxArray *prhs[4];   // Iteration Function (fcnName, niter/feval, fval/sse, x)
    // Statistics
    size_t ndec   = 0;
    // Settings
    bool enabled  = false;
} iterfun_cb_data_t;

// Termination Checks Structure
class OPTITermSettings
{
    public:
        OPTITermSettings(std::string solverNameIn) : solverName(solverNameIn) { }
        OPTITermSettings(std::string solverNameIn, double maxTimeIn, 
                                                   size_t maxIterIn,
                                                   size_t maxFevalIn,
                                                   size_t maxNodesIn) :
                                                   solverName(solverNameIn),
                                                   maxTime(maxTimeIn),
                                                   maxIter(maxIterIn),
                                                   maxFeval(maxFevalIn),
                                                   maxNodes(maxNodesIn) { }

        std::string solverName;
        double maxTime  = OPTI_MAX_TIME;
        size_t maxIter  = OPTI_MAX_ITER;
        size_t maxFeval = OPTI_MAX_FEVAL;
        size_t maxNodes = OPTI_MAX_NODES;
};

//
// Solver Licenses
//
enum class OptiSolverLicense {NONE, EPL, GPL, LGPL};

//
// License Information
//
#define OPTI_LGPL_LICENSE_LINK ("GNU Lesser General Public License: http://www.gnu.org/copyleft/lesser.html")
#define OPTI_GPL_LICENSE_LINK  ("GNU General Public License: https://www.gnu.org/copyleft/gpl.html")
#define OPTI_EPL_LICENSE_LINK  ("Eclipse Public License: http://opensource.org/licenses/eclipse-1.0")

//
// Solver Properties
//
class OptiSolverProperties
{
    public:
        OptiSolverProperties(std::string solverNameIn, std::string solverVersionIn, std::string solverLinkIn, OptiSolverLicense solverLicenseIn) :
                            solverName(solverNameIn),
                            solverVersion(solverVersionIn),
                            solverLink(solverLinkIn),
                            solverLicense(solverLicenseIn) { }
        OptiSolverProperties(std::string solverNameIn, double solverVersionIn, std::string solverLinkIn, OptiSolverLicense solverLicenseIn) :
                            solverName(solverNameIn),
                            solverLink(solverLinkIn),
                            solverLicense(solverLicenseIn) 
        { 
            char solverVersionStr[20];
            snprintf(solverVersionStr, 20, "%.2f", solverVersionIn);
            solverVersion = std::string(solverVersionStr);
        }
        std::string solverName;
        std::string solverVersion;
        std::string solverLink;
        OptiSolverLicense solverLicense;
};

//
// OPTI High Performance Timer
//
class OPTITimer
{
    public:
        OPTITimer(bool startNow = true);
        bool start(void);
        double execTime(void);

    private:
        LARGE_INTEGER _startCount;
        double _counterPeriod = 0.0;
        bool   _configured = false;
};

//
// Handy Dandy Class for Passing Around Mex Args + Solver Info
//
class OptiMexArgs
{
    public:
        OptiMexArgs(OptiSolverProperties solverInfoIn, const mxArray* prhsIn[], int nrhsIn, mxArray* plhsIn[], int nlhsIn) :
                    solverInfo(solverInfoIn),
                    prhs(prhsIn),
                    nrhs(nrhsIn),
                    plhs(plhsIn),
                    nlhs(nlhsIn) { }

        OptiSolverProperties solverInfo;
        const mxArray** prhs;
        int nrhs;
        mxArray** plhs;
        int nlhs;
};


//
// OPTI MEX Helper Methods
//
class OPTIMex
{
    public:
        // MATLAB Data Creation
        static mxArray* createDoubleScalar(double val = 0.0);
        static mxArray* createDoubleMatrix(size_t nrow, size_t ncol);
        static mxArray* createStruct(const char* fieldNames[], int numFields);
        static double* createFieldScalar(mxArray* data, const char* fieldName, double val = 0.0);
        static double* createFieldMatrix(mxArray* data, const char* fieldName, size_t nrow, size_t ncol);
        static mxArray* createFieldString(mxArray* data, const char* fieldName, const char* str);

        // Data Access
        static size_t getFieldNumElem(const mxArray* data, const char* fieldName);
        static size_t getNumElem(const mxArray* data);
        static size_t getNumRows(const mxArray* data);
        static size_t getNumCols(const mxArray* data);
        static mxArray* getField(const mxArray* data, const char* fieldName);
        static double*  getFieldPr(const mxArray* data, const char* fieldName);
        static std::string getFieldString(const mxArray* data, const char* fieldName);
        
        // Data Validation
        static bool isValidStruct(const mxArray* data);
        static bool isValidField(const mxArray* data, const char* field);
        static bool isDoubleScalar(const mxArray* data);
        static bool isDoubleVector(const mxArray* data);
        static bool isDoubleMatrix(const mxArray* data);
        static bool isString(const mxArray* data);

        static bool isRealDouble(const mxArray* data);
        static bool isScalar(const mxArray* data);
        static bool isVector(const mxArray* data);
        static bool isMatrix(const mxArray* data);        
        static bool isEmpty(const mxArray* data);

        // Options Structure Access
        static int getIntegerOption(const mxArray* opts, const char* optionName, int& option);
        static int getIntegerOption(const mxArray* opts, const char* optionName, size_t& option);
        static int getDoubleOption(const mxArray* opts, const char* optionName, double& option);

        // MATLAB Callback Data Init
        static void initMatlabCallbackData(matlab_cb_data_t& callbackData, const mxArray* problemData, size_t ndec, size_t ndata);
        static void initIterFunCallbackData(iterfun_cb_data_t& callbackData, mxArray* iterFunHandle, size_t ndec);

        // MATLAB Calling
        static double* callMatlabObjective(matlab_cb_data_t* mlData, double *x);
        static double* callMatlabGradient(matlab_cb_data_t* mlData, double *x);
        static bool    callMatlabIterFun(iterfun_cb_data_t* iterData, size_t niter, double fval, double* x);

        // Copy Helpers
        static void copyTranspose(double* dest, double* src, size_t nrows, size_t ncols);

        // Common Termination Checks
        static int checkTermination(OPTITimer& timer, size_t numIter, size_t numFeval, size_t numNodes, const OPTITermSettings& termSettings);
        
        // Error Reporting (throws exception too)
        static void error(const char* format, ...); 

        // Printing
        static void printSolverHeader(const OptiSolverProperties& solverInfo, const char* algorithmString, size_t ndec, size_t ndata);
        static void printSolverIter(size_t niter, double execTime, const char *header, const char *format, ...);
        static void printSolverInfo(const OptiSolverProperties& solver);

        // Ctrl-C Detection
        static bool ctrlCPressed(const char* solverName);

        // Standard Entry Check
        static bool checkNoArgs(const OptiMexArgs& args);
};

} // namespace opti_mex_utils



#endif  // OPTI_MEX_UTILS        