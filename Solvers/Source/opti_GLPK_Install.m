%% GLPK Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile GNU Linear Programming Kit (GLPK) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 1) Get GLPK
% GLPK is available from http://www.gnu.org/software/glpk/glpk.html. We 
% will create VS projects below.

% 2) Compile GLPK
% The easiest way to compile GLPK is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\glpk-4.48'; % FULL path to GLPK [MAX VER 4.48]
% opti_VSBuild('GLPK',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the GLPK MEX file. Once you have completed all the
% above steps, simply run this file to compile GLPK! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'glpkcc.cpp';
%Include Directories
inc = 'Include/Glpk';
%Lib Names [static libraries to link against]
libs = 'libglpk';
%Options
opts = [];
opts.verb = false;

%Compile
opti_solverMex('glpk',src,inc,libs,opts);