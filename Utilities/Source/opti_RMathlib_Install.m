%% R-Mathlub Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile the R Standalone Mathlib for use 
% with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 1) Get R source
% R is available from http://cran.stat.sfu.ca/. 
% Download the source.

% 2) Compile R standalone Math library
% The easiest way to compile R-Mathlib is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:
%
% % Visual Studio Builder Commands
% path = 'C:\Solvers\R-3.2.4';
% sdir = [path '\src\nmath']; 
% inc = [path '\src\include'];
% name = 'libRMathlib';
% opts = [];
% opts.compileAsCpp = true;
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','MATHLIB_STANDALONE'};
% opts.exPP = [opts.exPP 'HAVE_HYPOT','HAVE_EXPM1','HAVE_LOG1P']; %only if VS2015 or above!!
% opts.exclude = {'test.c'};
% VS_WriteProj(sdir,name,inc,opts)
% %%
% Once complete, you will have a directory called R\libRmathlib. Open the
% Visual Studio project file, then complete the following steps:
%   a) Copy from Utilities/Source/Include/RMathlib Rconfig.h and Rmath.h to
%   the R/src/nmath folder. 
%   b) Change nmath.h 110-112 to (VS2015)
%       - #define ML_POSINF	INFINITY
%       - #define ML_NEGINF	-INFINITY
%       - #define ML_NAN		NAN
%   OR VS2013
%       - #define ML_POSINF	std::numeric_limits<double>::infinity()
%       - #define ML_NEGINF	-std::numeric_limits<double>::infinity()
%       - #define ML_NAN std::numeric_limits<double>::quiet_NaN()
%   c) A number of files use <config.h> instead of "Rconfig.h". Compile and
%   look for the errors (pcauchy.c, fround.c, etc).
%   d) Comment line 70 in Arith.h (int R_finite(double);...)
%   e) VS2015 does not seem to like hexadecimal floating point constants.
%   Change line 86 in qbeta.c to DBL_1__eps    = 1 - DBL_EPSILON;
%   f) Comment line 80, 81 and 89 in Arith.h to remove macro redefinition
%   warning for ISNAN and R_FINITE
%   g) Add "#include "Rmath.h"" near the top of sunif.c (or ../Rmath.h)
%   h) Try compile, there should be a number of errors related to cannot
%   cast bool to Rboolean. For each error, manually cast the result (i.e.
%   add (Rboolean)... bit of a pain!).
%   i) Build a Win32 or x64 Release to compile the code.
%   j) Copy the generated .lib file to the following folder:
%
%   OPTI/Utilities/Source/lib/win32 or win64
%
%   You will also need to copy nmath.h to the following folder:
%
%   OPTI/Utilities/Source/Include/Rmathlib

% 3) RMathlib MEX Interface
% The R Math Library does not come with a MEX interface that I know of, so
% I wrote a simple one. Not all functions are present, but it works for my
% purposes. If there are functions you want added, let me know.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NLOPT MEX file. Once you have completed all 
% the above steps, simply run this file to compile NLOPT! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = 'rmathlibmex.cpp';
%Include Directories
inc = 'Include\Rmathlib';
%Lib Names [static libraries to link against]
libs = 'libRMathlib';
%Options
opts = [];
opts.verb = false;
opts.pp = 'MATHLIB_STANDALONE';
opts.util = true;

%Compile
opti_solverMex('rmathlib',src,inc,libs,opts);
