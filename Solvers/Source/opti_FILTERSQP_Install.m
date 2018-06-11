%% FILTERSQP Install for OPTI Toolbox
% Copyright (C) 2018 Jonathan Currie (IPL)

% This file will help you compile FILTER SQP for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2017
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% 1) Get FILTER SQP
% Unfortunately not publically available, but can be requested from the University of Dundee.

% 2) Compile FILTERSQP
% The easiest way to compile FILTERSQP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
path = 'C:\Solvers\filterSQP'; % FULL path to FilterSQP
opti_VSBuild('FilterSQP',path);