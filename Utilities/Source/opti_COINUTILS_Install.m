%% CoinUtils Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile CoinUtils for use with MATLAB.

% My build platform:
% - Windows 7 SP1 x64
% - Visual Studio 2013

% 1) Compile Clp as per opti_CLP_Install.m in Solvers/Source

% 2) Compile the MEX Files
% The code below will automatically include all required libraries and
% directories to build the MEX files. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!


%COIN READ
%MEX Interface Source Files
src = 'coinR.cpp';
%Include Directories
inc = {'..\..\Solvers\Source\Include\Coin','..\..\Solvers\Source\Include\Glpk'};
%Lib Names [static libraries to link against]
libs = {'libcoinutilsgmpl','libglpk'};
%Options
opts = [];
opts.verb = false;
opts.util = true;

%Compile
opti_solverMex('coinR',src,inc,libs,opts);

%COIN WRITE
%MEX Interface Source Files
src = 'coinW.cpp';
%Include Directories
inc = '..\..\Solvers\Source\Include\Coin';
%Lib Names [static libraries to link against]
libs = {'libcoinutils'};
%Options
opts = [];
opts.verb = false;
opts.util = true;

%Compile
opti_solverMex('coinW',src,inc,libs,opts);
