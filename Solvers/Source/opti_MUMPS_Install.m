%% MUMPS Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (IPL)

% This file will help you compile aMUltifrontal Massively Parallel sparse
% direct Solver (MUMPS) for use with MATLAB. 

% The supplied files and instructions are for compiling sequential double 
% precision MUMPS only.

% My build platform:
% - Windows 7 x64
% - Visual Studio 2013
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get MUMPS
% MUMPS is available from http://graal.ens-lyon.fr/MUMPS/. You will need to
% register before you can download.

% 2) Get METIS
% METIS is available from http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD.
% Download version 4.0.3 (last compatible version with MUMPS).

% 3) Compile MUMPS and METIS
% The easiest way to compile MUMPS is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\MUMPS_4.10.0'; % FULL path to MUMPS
% metispath = 'C:\Solvers\metis-4.0.3'; % FULL path to METIS [max version 4.0.3]
% opts = []; opts.expaths = metispath;
% opti_VSBuild('MUMPS',path,opts);

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MUMPS MEX file. Once you have completed all the
% above steps, simply run this file to compile MUMPS! You MUST BE in the 
% base directory of OPTI!

%Double MEX Interface Source Files
src = 'mumpsmex.c';
%Include Directories
inc = 'Include/Mumps';
%Lib Names [static libraries to link against]
libs = {'libdmumps_c','libdmumps_f','libpord','libseq_c','libseq_f','libmetis'};
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';
opts.pp = {'MUMPS_ARITH=2'};
opts.ifort = true;

%Compile
opti_solverMex('mumps',src,inc,libs,opts);

%Complex Double MEX Interface Source Files
src = 'mumpsmex.c';
%Include Directories
inc = 'Include/Mumps';
%Lib Names [static libraries to link against]
libs = {'libzmumps_c','libzmumps_f','libpord','libseq_c','libseq_f','libmetis'};
%Options
opts = [];
opts.verb = false;
opts.blas = 'MKL';
opts.pp = {'MUMPS_ARITH=8'};
opts.ifort = true;

%Compile
opti_solverMex('zmumpsmex',src,inc,libs,opts);

% METIS Reference:
% “A Fast and Highly Quality Multilevel Scheme for Partitioning Irregular 
% Graphs”. George Karypis and Vipin Kumar. SIAM Journal on Scientific 
% Computing, Vol. 20, No. 1, pp. 359—392, 1999.