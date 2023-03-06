%% OOQP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile Objective Orientated Quadratic
% Programming (OOQP) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get OOQP
% OOQP is available from http://pages.cs.wisc.edu/~swright/ooqp/. You will 
% need to register before you can download. We will create the VS project 
% below.

% 2) Compile OOQP
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Compilation Options (used when compiling the libs AND the mex file)
opts = [];
opts.pardiso = 'MKL'; %use Intel MKL's Pardiso Library (empty to not use pardiso)
opts.ma57 = 'Matlab'; %use MATLAB's supplied MA57 library (empty to not use ma57)

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\OOQP-0.99.22'; % FULL path to OOQP
% opti_VSBuild('OOQP',path,opts);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the OOQP MEX file. Once you have completed all the
% above steps, simply run this file to compile OOQP! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'ooqpmex.cpp';
%Include Directories
inc = 'Include/Ooqp';
%Lib Names [static libraries to link against]
libs = 'libooqp';
%Options (note options from above used here too)
opts.verb = false;
opts.blas = 'mkl'; %note only compatible with MKL as BLAS calls are of the form dscal_

%Compile
opti_solverMex('ooqp',src,inc,libs,opts);