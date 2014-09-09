%% MUMPS Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile aMUltifrontal Massively Parallel sparse
% direct Solver (MUMPS) for use with MATLAB. 

% The supplied files and instructions are for compiling sequential double 
% precision MUMPS only.

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
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

mumpspath = 'C:\Solvers\MUMPS_4.10.0'; % FULL path to MUMPS
metispath = 'C:\Solvers\metis-4.0.3'; % FULL path to METIS

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('MUMPS',{mumpspath,metispath},cd);

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MUMPS MEX file. Once you have completed all the
% above steps, simply run this file to compile MUMPS! You MUST BE in the 
% base directory of OPTI!

clear mumps

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('MUMPS MEX FILE INSTALL\n\n');

%Get MUMPS Libraries
post = [' -IInclude/Mumps -L' libdir ' -llibdmumps_c -llibdmumps_f -llibpord -llibseq_c -llibseq_f -llibmetis -DMUMPS_ARITH=2'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common
post = [post ' -output mumps'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims mumpsmex.c';
try
    eval([pre post])
    movefile(['mumps.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:mumps','Error Compiling MUMPS!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');

% METIS Reference:
% “A Fast and Highly Quality Multilevel Scheme for Partitioning Irregular 
% Graphs”. George Karypis and Vipin Kumar. SIAM Journal on Scientific 
% Computing, Vol. 20, No. 1, pp. 359—392, 1999.

% Optional ZMUMPS Install
clear zmumpsmex

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('ZMUMPS MEX FILE INSTALL\n\n');

%Get MUMPS Libraries
post = [' -IInclude/Mumps -L' libdir ' -llibzmumps_c -llibzmumps_f -llibpord -llibseq_c -llibseq_f -llibmetis -DMUMPS_ARITH=8'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common
post = [post ' -output zmumpsmex'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims mumpsmex.c';
try
    eval([pre post])
    movefile(['zmumpsmex.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:mumps','Error Compiling ZMUMPS!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');