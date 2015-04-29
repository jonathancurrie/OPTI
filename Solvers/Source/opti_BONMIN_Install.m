%% BONMIN Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Basic Open-source Nonlinear Mixed INteger
% programming (BONMIN) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 0) Complete Compilation as per OPTI instructions for CLP and CBC.

% 1) Get BONMIN
% BONMIN is available from http://www.coin-or.org/Bonmin/. Download 
% the source. 

% 2) Compile BONMIN
% The easiest way to compile BONMIN is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Compilation Options (used when compiling the libs AND the mex file) - see opti_VSBuild for details
opts = [];
opts.pardiso = ''; %do not include PARDISO (too big for general distribution)
opts.ma57 = 'Matlab'; %use MATLAB's supplied MA57 library (empty to not use ma57)
opts.mumps = true; %link MUMPS (add path below when compiling lib)
opts.ma27 = ''; %do not link ma27

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\Bonmin-1.8.2\Bonmin'; %FULL path to BONMIN
% ipoptpath = 'C:\Solvers\Ipopt-3.12.3\Ipopt'; %FULL path to IPOPT
% cbcpath = 'C:\Solvers\Cbc-2.9.4\Cbc'; % FULL path to CBC
% clppath = 'C:\Solvers\Clp-1.16.6\Clp'; % FULL path to CLP
% mumpspath = 'C:\Solvers\MUMPS_4.10.0'; %FULL path to MUMPS (or leave blank to skip linking MUMPS)
% metispath = 'C:\Solvers\metis-4.0.3'; % FULL path to METIS (leave blank if not linking MUMPS) [max version 4.0.3]
% opts.expaths = {ipoptpath,cbcpath,clppath,mumpspath,metispath};
% opti_VSBuild('BONMIN',path,opts);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the BONMIN MEX file. Once you have completed all 
% the above steps, simply run this file to compile BONMIN! You MUST BE in 
% the base directory of OPTI!

%MEX Interface Source Files
src = {'bonmin/matlabexception.cpp','bonmin/matlabfunctionhandle.cpp','bonmin/matlabjournal.cpp'...
       'bonmin/iterate.cpp','bonmin/bonminoptions.cpp','bonmin/options.cpp','bonmin/sparsematrix.cpp',...
       'bonmin/callbackfunctions.cpp','bonmin/matlabinfo.cpp','bonmin/matlabprogram.cpp','bonmin/bonminmex.cpp'};
%Include Directories
inc = {'Include/bonmin','bonmin/Include','Include/Ipopt','ipopt/Include','Include/BuildTools','Include/Cbc',...
        'Include/Osi','Include/Cgl','Include/Clp','Include/Coin'};
%Lib Names [static libraries to link against]
libs = {'libbonmin','libipoptbm','libcbc','libosi','libcgl','libclp','libcoinutils'};
%Options (note options from above used here too)
opts.verb = false;
opts.pp = {'BONMIN_BUILD','IPOPT_BUILD'};
opts.blas = 'MKL'; 

%Compile
opti_solverMex('bonmin',src,inc,libs,opts);
