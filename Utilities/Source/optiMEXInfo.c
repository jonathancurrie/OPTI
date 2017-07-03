/*
 * OPTI MEX Info - Returns information about the build information of the OPTI MEX Files
 */
#include "mex.h"
#include "opti_util.h"

#define STR_EXPAND(tok) #tok
#define str(tok) STR_EXPAND(tok)

void printBuildInfo(void);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs == 0)
    {
        printBuildInfo();
        return;
    }
    else
    {
        plhs[0] = mxCreateString(str(MEX_HASH));
        plhs[1] = mxCreateString(__DATE__);
        plhs[2] = mxCreateString(__TIME__);
    }
}


void printBuildInfo(void)
{
    char vbuf[10]; getVSVer(vbuf);  
    mexPrintf("\n-----------------------------------------------------------\n");
    mexPrintf(" OPTI MEX Info [Built %s, %s, VS%s, MATLAB %s]\n",__TIME__,__DATE__,vbuf, str(ML_VER));
    mexPrintf(" SHA-256: %s\n", str(MEX_HASH));
    mexPrintf("-----------------------------------------------------------\n");
}