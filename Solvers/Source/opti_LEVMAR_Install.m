%% LEVMAR Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (Control Engineering)

% This file will help you compile Levenberg-Marquardt in C/C++ (LEVMAR) for 
% use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get LEVMAR
% LEVMAR is available from http://www.ics.forth.gr/~lourakis/levmar/. We 
% will create the VS project below.

% 2) Compile LEVMAR
% The easiest way to compile LEVMAR is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

% Intel MKL 11.3 R3 Change (due to workspace query NULLs) in Axb_core.c line 1026:
% if(!nb){
%     LM_REAL tmp;
% 	  LM_REAL tmpa;
% 	  LM_REAL tmppiv;
% 
%     work_sz=-1; // workspace query; optimal size is returned in tmp
%     SYTRF("L", (int *)&m, &tmpa, (int *)&m, &tmppiv, (LM_REAL *)&tmp, (int *)&work_sz, (int *)&info);
%     nb=((int)tmp)/m; // optimal worksize is m*nb
%   }

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\levmar-2.6'; %FULL path to LEVMAR
% opti_VSBuild('LEVMAR',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LEVMAR MEX file. Once you have completed all the
% above steps, simply run this file to compile LEVMAR! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'levmarmex.c';
%Include Directories
inc = {'Include\Levmar'};
%Lib Names [static libraries to link against]
libs = 'liblevmar';
%Options
opts = [];
opts.verb = false;
opts.blas = 'mkl';

%Compile
opti_solverMex('levmar',src,inc,libs,opts);
