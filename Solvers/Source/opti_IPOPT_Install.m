%% IPOPT Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile Interior Point OPTimizer (IPOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2015
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% NOTE - From OPTI v1.71 IPOPT can now be dynamically linked against the
% MathWorks supplied libmwma57.dll (HSL MA57). This step is optional
% as we are also compiling MUMPS. Alternatively you can skip MUMPS, 
% and just use MA57! Also be aware libmwma57.dll does not play well
% on unconstrained problems due in part to missing MeTiS, thus the 
% ma57 pivot order option is overidden automatically.

% To recompile you will need to get / do the following:

% 1) Get IPOPT
% IPOPT is available from http://www.coin-or.org/download/source/Ipopt/.
% Download the latest version of the source.

% 2) Compile IPOPT
% The easiest way to compile IPOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Compilation Options (used when compiling the libs AND the mex file) - see opti_VSBuild for details
opts = [];
opts.pardiso = 'MKL'; %use Intel MKL's Pardiso Library (empty to not use pardiso)
opts.ma57 = 'Matlab'; %use MATLAB's supplied MA57 library (empty to not use ma57)
opts.mumps = true; %link MUMPS (add path below when compiling lib)
opts.ma27 = ''; %do not link ma27
opts.linloader = false; %do not use HSL's linear solver dynamically loaded library (you must compile this separately)

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\Ipopt-3.12.6\Ipopt'; %FULL path to IPOPT
% mumpspath = 'C:\Solvers\MUMPS_4.10.0'; %FULL path to MUMPS (or leave blank to skip linking MUMPS)
% metispath = 'C:\Solvers\metis-4.0.3'; % FULL path to METIS (leave blank if not linking MUMPS) [max version 4.0.3]
% opts.expaths = {mumpspath,metispath};
% opti_VSBuild('IPOPT',path,opts);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the IPOPT MEX file. Once you have completed all the
% above steps, simply run this file to compile IPOPT! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = {'ipopt/matlabexception.cpp','ipopt/matlabfunctionhandle.cpp','ipopt/matlabjournal.cpp'...
       'ipopt/iterate.cpp','ipopt/ipoptoptions.cpp','ipopt/options.cpp','ipopt/sparsematrix.cpp',...
       'ipopt/callbackfunctions.cpp','ipopt/matlabinfo.cpp','ipopt/matlabprogram.cpp','ipopt/ipopt.cpp'};
%Include Directories
inc = {'Include/Ipopt','ipopt/Include','Include/BuildTools'};
%Lib Names [static libraries to link against]
libs = 'libipopt';
%Options (note options from above used here too)
opts.verb = false;
opts.pp = 'IPOPT_BUILD';
opts.blas = 'MKL'; 

%Compile
opti_solverMex('ipopt',src,inc,libs,opts);

%% Reproducible Results with IPOPT
% By default IPOPT is compiled above with the multi-threaded MKL libraries.
% Due to threading order affecting the order of floating point operations,
% consecutive runs of the SAME problem may produce DIFFERENT results. This
% effect is problem and platform dependent.

% To obtain reproducible results, you must remove all multi-threading from
% IPOPT. The easiest way to do this is to only use MUMPS, and compile
% against MKL_SEQ or NETLIB as the BLAS library. Do not use MA57 or
% PARDISO. Note however you will have a substantial speed-hit for doing
% this.





