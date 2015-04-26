%% PSWARM Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile PSwarm for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get PSwarm
% PSwarm is available from http://www.norg.uminho.pt/aivaz/pswarm/. 
% Download PPSwarm_vxx.zip (C Version) and unzip to a suitable location.

% 2) Compile PSwarm
% The easiest way to compile PSwarm is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\PPSwarm_v1_5'; % FULL path to PSwarm
% opti_VSBuild('PSwarm',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the PSwarm MEX file. Once you have completed all the
% above steps, simply run this file to compile PSwarm! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'pswarmmex.c';
%Include Directories
inc = 'Include/Pswarm';
%Lib Names [static libraries to link against]
libs = 'libpswarm';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';

%Compile
opti_solverMex('pswarm',src,inc,libs,opts);