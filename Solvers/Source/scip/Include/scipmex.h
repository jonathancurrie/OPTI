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

//Error catching macro
#define SCIP_ERR(x,msg) if(x != SCIP_OKAY) {sprintf(msgbuf,"%s, Error Code: %d",msg,x); mexErrMsgTxt(msgbuf);}

//Add Ctrl-C event handler
SCIP_RETCODE SCIPincludeCtrlCEventHdlr(SCIP* scip);

//Add Nonlinear constraint / objective
double addNonlinearCon(SCIP* scip, SCIP_VAR** vars, double *instr, size_t no_instr, double lhs, double rhs, double *x0, int nlno, bool isObj);

#endif