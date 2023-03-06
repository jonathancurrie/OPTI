%% GSL Install for OPTI Toolbox
% Copyright (C) 2017 Jonathan Currie (Control Engineering)
clc
% This file will help you compile GNU Scientific Library for use with 
% MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2017

% To recompile you will need to get / do the following:

% 1) Get GSL
% GSL for Windows is available from https://github.com/BrianGladman/gsl

% 2) Compile GSL
% The easiest way to compile GSL is to use the Visual Studio Projects in
% the above download. Open the gsl.lib solution, build the "gslhdrs"
% project, the build the "gsllib" project. Copy the resulting .lib and
% header files (/lib/x64/Release and /gsl).

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the GSL MEX file. Once you have completed all 
% the above steps, simply run this file to compile GSL! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = {'gsl/gslmex.cpp', 'gsl/gslmex_nls.cpp'};
%Include Directories
inc = [];
%Lib Names [static libraries to link against]
libs = {'libgsl'};
%Options
opts = [];
opts.pp = {};
opts.blas = 'MKL'; 
opts.verb = false;

%Compile
opti_solverMex('gsl',src,inc,libs,opts);
