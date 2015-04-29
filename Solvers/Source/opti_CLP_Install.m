%% CLP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)
clc
% This file will help you compile Coin-Or Linear Programming for use with 
% MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% To recompile you will need to get / do the following:

% 1) Get CLP
% CLP is available from http://www.coin-or.org/projects/Clp.xml. Download 
% the source. Unzip to a suitable temporary location.

% 2) Compile CLP & COIN Utils
% The easiest way to compile CLP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required paths on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\Clp-1.16.6\Clp'; % FULL path to CLP
% glpkpath = 'C:\Solvers\glpk-4.48'; % FULL path to GLPK (or leave blank [MAX VER 4.48])
% opts = []; opts.expath = glpkpath;
% opti_VSBuild('CLP',path,opts);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CLP MEX file. Once you have completed all 
% the above steps, simply run this file to compile CLP! You MUST BE in 
% the base directory of OPTI!

%Enable for Aboca Build (must compile using Intel C++ & VS2012/2013 Linker)
haveABC = false;

%MEX Interface Source Files
src = 'clpmex.cpp';
%Include Directories
inc = {'Include/Clp','Include/Coin'};
%Lib Names [static libraries to link against]
libs = {'libCoinUtils','libosi'};
%Options
opts = [];
opts.pp = {'COIN_MSVS'};
opts.verb = false;

%Optional Aboca Setup
if(haveABC)
    opts.pp = [opts.pp, 'INTEL_COMPILER', 'CLP_HAS_ABC=4', '__BYTE_ORDER=__LITTLE_ENDIAN'];
    libs = [libs, 'libclpabc'];
else
    libs = [libs, 'libclp'];
end

%Compile
opti_solverMex('clp',src,inc,libs,opts);
