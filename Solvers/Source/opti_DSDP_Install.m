%% DSDP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile DSDP for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get DSDP
% DSDP is available from http://www.mcs.anl.gov/hs/software/DSDP/. Download
% it and unzip the directory.

% 2) Compile DSDP
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\DSDP5.8'; % FULL path to DSDP
% opti_VSBuild('DSDP',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the DSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile DSDP! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'dsdpmex.c';
%Include Directories
inc = 'Include/DSDP';
%Lib Names [static libraries to link against]
libs = 'libdsdp';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';

%Compile
opti_solverMex('dsdp',src,inc,libs,opts);