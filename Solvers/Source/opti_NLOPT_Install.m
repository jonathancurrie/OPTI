%% NLOPT Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile NonLinear OPTimization (NLOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 SP1 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 1) Get NLOPT
% NLOPT is available from http://ab-initio.mit.edu/wiki/index.php/NLopt. 
% Download the source.

% 2) Compile NLOPT
% The easiest way to compile NLOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\nlopt-2.4.2'; % FULL path to NLOPT
% opti_VSBuild('NLOPT',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NLOPT MEX file. Once you have completed all 
% the above steps, simply run this file to compile NLOPT! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = 'nloptmex.c';
%Include Directories
inc = {'Include\Nlopt','nlopt\Include'};
%Lib Names [static libraries to link against]
libs = 'libnlopt';
%Options
opts = [];
opts.verb = false;

%Compile
opti_solverMex('nlopt',src,inc,libs,opts);

