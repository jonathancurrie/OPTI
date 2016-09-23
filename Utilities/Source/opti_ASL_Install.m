%% ASL Install for OPTI Toolbox
% Supplied binaries are built from Netlib's AMPL Solver Library Interface

%   Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile AMPL Solver Library (ASL) for use with 
% MATLAB.

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013

% 1) Get AMPL Solver Library
% The generic NL reader for AMPL is available free from AMPL 
% (http://www.ampl.com/netlib/ampl/solvers/). Download the gzipped tar of
% all files AND degree.c individually if it is empty! Note this is not the 
% AMPL engine which is a commerical product, but code to allow people to 
% connect their solvers to AMPL.

% 2) Compile ASL
% The easiest way to compile ASL is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\ASL'; % FULL path to ASL
% opti_VSBuild('ASL',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MEX file. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = 'amplmex.c';
%Include Directories
inc = 'Include/Asl';
%Lib Names [static libraries to link against]
libs = 'libasl';
%Options
opts = [];
opts.verb = false;
opts.debug = false;
opts.util = true;
opts.pp = 'NO_STDIO1';

%Compile
opti_solverMex('asl',src,inc,libs,opts);