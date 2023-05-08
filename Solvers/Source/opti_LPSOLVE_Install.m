%% LP_SOLVE Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (Control Engineering)

% This file will help you compile LP_SOLVE for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 1) Get LP_SOLVE
% LP_SOLVE is available from
% http://sourceforge.net/projects/lpsolve/files/lpsolve/. Download the
% source (lp_solve_5.5.2.0_source.tar.gz) or later version.

% 2) Compile LP_SOLVE
% The easiest way to compile LP_SOLVE is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\lp_solve_5.5'; % FULL path to LP_SOLVE
% opti_VSBuild('LPSOLVE',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LP_SOLVE MEX file. Once you have completed all 
% the above steps, simply run this file to compile LP_SOLVE! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = {'lpsolve/lpsolve.c','lpsolve/matlab.c'};
%Include Directories
inc = {'Include/Lpsolve','lpsolve/Include'};
%Lib Names [static libraries to link against]
libs = 'liblpsolve';
%Options
opts = [];
opts.verb = false;
opts.pp = {'LPSOLVEAPIFROMLIB','MATLAB','WIN32'};

%Compile
opti_solverMex('lp_solve',src,inc,libs,opts);