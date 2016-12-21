%% MINPACK Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile MINPACK HYBRJ + HYBRJ + LMDER + LMDIF for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get MINPACK
% HYBRJ + LMDER are part of MINPACK available from 
% http://www.netlib.org/minpack/. Download and combine lmder, lmdir, hybrd 
% + hybrj + all dependencies into single directory.

% 2) Compile HYBRJ + HYBRD + LMDER + LMDIF
% The easiest way to compile MINPACK is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\minpack'; % FULL path to MINPACK
% opti_VSBuild('minpack',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the HYBRJ MEX file. Once you have completed all the
% above steps, simply run this file to compile HYBRJ! You MUST BE in the 
% base directory of OPTI!

%HBYRJ MEX Interface Source Files
src = 'hybrjmex.c';
%Lib Names [static libraries to link against]
libs = 'libminpack';
%Options
opts = [];
opts.verb = false;
opts.ifort = true;

%Compile
opti_solverMex('hybrj',src,[],libs,opts);

%LMDER MEX Interface Source Files
src = 'lmdermex.c';
%Lib Names [static libraries to link against]
libs = 'libminpack';
%Options
opts = [];
opts.verb = false;
opts.ifort = true;

%Compile
opti_solverMex('lmder',src,[],libs,opts);