/* SCIPMEX - A MATLAB MEX Interface to GSL
 * Released Under the BSD 3-Clause License:
 * https://www.inverseproblem.co.nz/OPTI/index.php/DL/License
 *
 * Copyright (C) Jonathan Currie 2017
 * www.inverseproblem.co.nz
 */

#pragma once

#include "mex.h"
#include "gsl_version.h"

// Problem Type Enum
enum class GslProbType {UNKNOWN,NLS};

//Argument Enums (in expected order of arguments)
enum {ePROB, eOPTS};                   

//PRHS Defines    
#define pPROB   prhs[ePROB]
#define pOPTS   prhs[eOPTS]

// Solve a NLS using GSL
void gslSolveNLS(const mxArray *prhs[], int nrhs, mxArray *plhs[], int nlhs);