%% FILTERSD Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile FILTER SD for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% 1) Get FILTER SD
% FilterSD is available from http://www.coin-or.org/download/source/filterSD/
% Download the .zip file for use on Windows.

% 2) Compile FILTERSD + FILTERSDSP
% The easiest way to compile FILTERSD is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\filterSD-1.0.0'; % FULL path to FilterSD [note I have modified quite a bit of the filterSD source for OPTI]
% opti_VSBuild('FilterSD',path);

% 3) FILTERSD MEX Interface
% The FILTERSD MEX Interface is a simple MEX interface I wrote to use
% FILTERSD.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the FILTERSD MEX file. Once you have completed all the
% above steps, simply run this file to compile FILTERSD! You MUST BE in the 
% base directory of OPTI!

%FILTERSD Dense MEX Interface Source Files
src = 'filtersdmex.c';
%Lib Names [static libraries to link against]
libs = 'libfiltersd';
%Options
opts = [];
opts.verb = false;
opts.ifort = true;

%Compile
opti_solverMex('filtersd',src,[],libs,opts);


%FILTERSD Sparse MEX Interface Source Files
src = 'filtersdmex.c';
%Lib Names [static libraries to link against]
libs = 'libfiltersdsp';
%Options
opts = [];
opts.verb = false;
opts.pp = {'SPARSEVER'};
opts.ifort = true;

%Compile
opti_solverMex('filtersdsp',src,[],libs,opts);