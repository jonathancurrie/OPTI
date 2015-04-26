/* SCIPMEX - A MATLAB MEX Interface to SCIP
 * Released Under the BSD 3-Clause License:
 * http://www.i2c2.aut.ac.nz/Wiki/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2013
 * www.i2c2.aut.ac.nz
 */

#ifndef SCIPMEXINC
#define SCIPMEXINC

#include "scip/scip.h"
#include "mex.h"

//Message buffer size
#define BUFSIZE 1024
//Function buffer size [maximum expressions & variables to hold for post-processing]
#define MAX_DEPTH 512
//Enable for debug print out
// #define DEBUG 1

//Convert Error Code to Message (from http://scip.zib.de/doc-3.0.2/html_devel/type__retcode_8h_source.shtml)
char *scipErrCode(int x)
{
    switch(x)
    {
        case SCIP_OKAY: return "Normal Termination";
        case SCIP_ERROR: return "Unspecified Error";
        case SCIP_NOMEMORY: return "Insufficient Memory Error";
        case SCIP_READERROR: return "Read Error";
        case SCIP_WRITEERROR: return "Write Error";
        case SCIP_NOFILE: return "File Not Found Error";
        case SCIP_FILECREATEERROR: return "Cannot Create File";
        case SCIP_LPERROR: return "Error in LP Solver";
        case SCIP_NOPROBLEM: return "No Problem Exists";
        case SCIP_INVALIDCALL: return "Method Cannot Be Called at This Time in Solution Process";
        case SCIP_INVALIDDATA: return "Error In Input Data";
        case SCIP_INVALIDRESULT: return "Method Returned An Invalid Result Code";
        case SCIP_PLUGINNOTFOUND: return "A required plugin was not found";
        case SCIP_PARAMETERUNKNOWN: return "The parameter with the given name was not found";
        case SCIP_PARAMETERWRONGTYPE: return "The parameter is not of the expected type";
        case SCIP_PARAMETERWRONGVAL: return "The value is invalid for the given parameter";
        case SCIP_KEYALREADYEXISTING: return "The given key is already existing in table";
        case SCIP_MAXDEPTHLEVEL: return "Maximal branching depth level exceeded";
        case SCIP_BRANCHERROR: return "No branching could be created";
        default: return "Unknown Error Code";
    }
}

//Error catching macro
#define SCIP_ERR(x,msg) if(x != SCIP_OKAY) {sprintf(msgbuf,"%s, Error %s (Code: %d)",msg,scipErrCode(x),x); mexErrMsgTxt(msgbuf);}

//Add Ctrl-C event handler
SCIP_RETCODE SCIPincludeCtrlCEventHdlr(SCIP* scip);

//Add Nonlinear constraint / objective
double addNonlinearCon(SCIP* scip, SCIP_VAR** vars, double *instr, size_t no_instr, double lhs, double rhs, double *x0, int nlno, bool isObj);

#endif