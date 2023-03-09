%% L-BFGS-B Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile Limited Memory Broyden-Fletcher-Goldfarb-
% Shanno Bounded Optimization (L-BFGS-B) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get L-BFGS-B
% L-BFGS-B is available from 
% http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html. We will create 
% VS projects below.

% 2) Compile L-BFGS-B
% The easiest way to compile L-BFGS-B is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\Lbfgsb.3.0'; % FULL path to L-BFGS-B
% opti_VSBuild('lbfgsb',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the L-BFGS-B MEX file. Once you have completed all the
% above steps, simply run this file to compile L-BFGS-B! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = {'lbfgsb/lbfgsb.cpp','lbfgsb/lbfgsb_program.cpp','lbfgsb/program.cpp'};
%Include Directories
inc = 'lbfgsb/Include';
%Lib Names [static libraries to link against]
libs = 'liblbfgsb';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';
opts.ifort = true;

%Compile
opti_solverMex('lbfgsb',src,inc,libs,opts);