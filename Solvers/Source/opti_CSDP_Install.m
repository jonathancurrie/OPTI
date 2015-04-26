%% CSDP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile CSDP for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get CSDP
% CSDP is available from https://projects.coin-or.org/Csdp/. Download
% it and unzip the directory.

% 2) Compile CSDP
% The easiest way to compile CSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\CSDP 6.2beta'; % FULL path to CSDP
% opti_VSBuild('CSDP',path); %note I have modified some of the CSDP source files, see opti_VSBuild for details

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile CSDP! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'csdpmex.c';
%Include Directories
inc = 'Include/Csdp';
%Lib Names [static libraries to link against]
libs = 'libcsdp';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';
opts.pp = {'NOSHORTS'};
opts.expre = 'LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:vcompd.lib /NODEFAULTLIB:vcomp.lib"'; %don't link against default VC++ OpenMP Lib, use Intel One (remove if not using Intel MKL)

%Compile
opti_solverMex('csdp',src,inc,libs,opts);