%% MKL Trust Region NLS Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile the Intel Math Kernel Library (MKL) 
% Trust Region NLS function for use with MATLAB. NOTE you must NOT link the 
% threaded MKL libraries as the MATLAB callback function is not thread safe!

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get and Install Intel MKL
% http://software.intel.com/en-us/articles/intel-mkl/

% 2) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MKL TR NLS MEX file. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'mkltrnls.c';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL_SEQ'; %note sequential MKL!

%Compile
opti_solverMex('mkltrnls',src,[],[],opts);