%% NL2SOL Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile NL2SOL for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get NL2SOL
% NL2SOL is available in multiple variants, however the most recent is 
% included in the PORT library, available from:
% http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/port
% The library contains functions for a range of mathematical functions,
% however we will just be using the NL2SOL variants.

% 2) Compile NL2SOL + NL2SNO (DN2F, DN2G, DN2FB, DN2GB)
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\port'; % FULL path to NL2SOL
% opti_VSBuild('NL2SOL',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NL2SOL MEX file. Once you have completed all the
% above steps, simply run this file to compile NL2SOL! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'nl2solmex.c';
%Lib Names [static libraries to link against]
libs = 'libnl2sol';
%Options
opts = [];
opts.verb = false;
opts.ifort = true; %requires Intel Fortran Dynamic Libraries

%Compile
opti_solverMex('nl2sol',src,[],libs,opts);