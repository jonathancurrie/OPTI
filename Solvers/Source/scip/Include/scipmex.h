/* SCIPMEX - A MATLAB MEX Interface to SCIP
 * Released Under the BSD 3-Clause License:
 * https://www.controlengineering.co.nz/Wikis/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2018
 * www.inverseproblem.co.nz
 */

#ifndef SCIPMEXINC
#define SCIPMEXINC

#include "scip/scip.h"
#include "mex.h"

//Message buffer size
#define BUFSIZE 2048
//Function buffer size [maximum expressions & variables to hold for post-processing]
#define MAX_DEPTH 512
//Enable for debug print out
// #define DEBUG 1

//Convert Error Code to Message (from http://scip.zib.de/doc-3.0.2/html_devel/type__retcode_8h_source.shtml)
char *scipErrCode(int x);

//Error catching macro
#define SCIP_ERR(rc,msg) if(rc != SCIP_OKAY) {sprintf(msgbuf,"%s, Error Code: %d",msg,rc); mexErrMsgTxt(msgbuf);}

//Add Ctrl-C event handler
SCIP_RETCODE SCIPincludeCtrlCEventHdlr(SCIP* scip);

//Add Nonlinear constraint / objective
double addNonlinearCon(SCIP* scip, SCIP_VAR** vars, double *instr, size_t no_instr, double lhs, double rhs, double *x0, size_t nlno, bool isObj);

#endif