%% M1QN3 Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile M1QN3 for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get M1QN3
% M1QN3 is available from 
% https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html
% Download and Unzip the folder, we will create a Visual Studio project below.

% 2) Compile M1QN3
% The easiest way to compile M1QN3 is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\m1qn3-3.3-distrib'; % FULL path to M1QN3
% opti_VSBuild('M1QN3',path);

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the M1QN3 MEX file. Once you have completed all the
% above steps, simply run this file to compile M1QN3! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'm1qn3mex.c';
%Lib Names [static libraries to link against]
libs = 'libm1qn3';
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';
opts.ifort = true;

%Compile
opti_solverMex('m1qn3',src,[],libs,opts);