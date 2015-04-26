%% CBC Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Coin-Or Branch and Cut for use with 
% MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 0) Compile CLP as per the instructions in opti_CLP_Install.m

% 1) Get CBC
% CBC is available from https://projects.coin-or.org/Cbc. Download 
% the source.

% 2) Compile CBC
% The easiest way to compile CBC is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\Cbc-2.9.3\Cbc'; % FULL path to CBC
% opti_VSBuild('CBC',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CBC MEX file. Once you have completed all 
% the above steps, simply run this file to compile CBC! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = 'cbcmex.cpp';
%Include Directories
inc = {'Include/Cbc','Include/Osi','Include/Cgl','Include/Clp','Include/Coin'};
%Lib Names [static libraries to link against]
libs = {'libcbc','libosi','libcgl','libclp','libcoinutils'};
%Options
opts = [];
opts.verb = false;
opts.pp = {'COIN_MSVS'};

%Compile
opti_solverMex('cbc',src,inc,libs,opts);