// OPTI Utilities Header File
// Collection of functions and macros for helping with OPTI MEX files

#ifndef OPTI_UTILS_HDR
#define OPTI_UTILS_HDR

#if (_MSC_VER == 1500)
    #define VS_VER "2008"
#elif (_MSC_VER == 1600)
    #define VS_VER "2010"
#elif (_MSC_VER == 1700)
    #define VS_VER "2012"
#elif (_MSC_VER == 1800)
    #define VS_VER "2013"
#elif (_MSC_VER == 1900)
    #define VS_VER "2015"
#elif (_MSC_VER == 1910)
    #define VS_VER "2017"
#else
    #define VS_VER "?"
#endif

#ifndef ML_VER
#error Please define the MATLAB version when building the MEX file
#endif
#ifndef OPTI_VER
#error Please define the OPTI version when building the MEX file
#endif

#define STR_EXPAND_OPTI(tok) #tok
#define stringify(tok) STR_EXPAND_OPTI(tok)
#define PRINT_BUILD_INFO mexPrintf("  - Built %s, Visual Studio %s, MATLAB %s, OPTI v%s\n", __DATE__, VS_VER, stringify(ML_VER), stringify(OPTI_VER));


#include "mex.h"
#include <string.h>

inline void CheckOptiVersion(const mxArray *opts)
{
    static bool displayedWarning = false;
    double localVer = 0.0;
    if (opts != NULL)
    {
        if(mxIsStruct(opts) && mxGetField(opts,0,"optiver") && !mxIsEmpty(mxGetField(opts,0,"optiver")))
        {
            localVer = *mxGetPr(mxGetField(opts,0,"optiver"));
            if(OPTI_VER != localVer)
            {
                if (displayedWarning == false)
                {
                    char buf[256];
                    sprintf(buf, "The MEX File Version (%.2f) does not match OPTI's Version (%.2f), please run opti_Install.m to update your MEX files.", OPTI_VER, localVer); 
                    mexWarnMsgTxt(buf);
                    displayedWarning = true; // don't spam
                }
            }
        }
    }
}

// OPTI Helper Methods
inline void getIntegerOption(const mxArray *opts, const char *name, int *var)
{
    if (opts == NULL) 
    { 
        mexErrMsgTxt("Supplied Options mxArray was NULL!\n");
        return; 
    }
    
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
    {
        *var = (int)*mxGetPr(mxGetField(opts,0,name));
    }
}
inline void getDoubleOption(const mxArray *opts, const char *name, double *var)
{
    if (opts == NULL) 
    { 
        mexErrMsgTxt("Supplied Options mxArray was NULL!\n");
        return; 
    }
    
    if(mxGetField(opts,0,name) && !mxIsEmpty(mxGetField(opts,0,name)))
    {
        *var = *mxGetPr(mxGetField(opts,0,name));
    }
}


// Error Message with Printf Functionality
#define ERR_BUF_LEN (1024)
#define MEX_ERR(fmt, ...) \
do { \
    char errBuf[ERR_BUF_LEN]; \
   snprintf(errBuf, ERR_BUF_LEN, fmt, __VA_ARGS__); \
   mexErrMsgTxt(errBuf); \
} while(0); 
        

// Ctrl-C Detection 
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
        
#endif  // OPTI_UTILS_HDR        