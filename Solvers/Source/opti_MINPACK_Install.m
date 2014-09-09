%% MINPACK Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile MINPACK HYBRJ + HYBRJ + LMDER + LMDIF for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get MINPACK
% HYBRJ + LMDER are part of MINPACK available from 
% http://www.netlib.org/minpack/. Download and combine lmder, lmdir, hybrd 
% + hybrj + all dependencies into single directory.

% 2) Compile HYBRJ + HYBRD + LMDER + LMDIF
% The easiest way to compile MINPACK is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

minpkpath = 'C:\Solvers\minpack'; % FULL path to MINPACK

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('MINPACK',minpkpath,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the HYBRJ MEX file. Once you have completed all the
% above steps, simply run this file to compile HYBRJ! You MUST BE in the 
% base directory of OPTI!

clear hybrj lmder

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('HYBRJ + HYBRD MEX FILE INSTALL\n\n');

%Get HYBRJ + HYBRD Libraries
post = [' -L' libdir ' -llibminpack -llibut -output hybrj'];
%Get Intel Fortran Libraries (for HYBRJ build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims hybrjmex.c';
try
    eval([pre post])
    movefile(['hybrj.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lmder','Error Compiling HYBRJ / HYBRD!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');

fprintf('\n------------------------------------------------\n');
fprintf('LM_DER + LM_DIF MEX FILE INSTALL\n\n');

%Get LMDER + LMDIF Libraries
post = [' -L' libdir ' -llibminpack -llibut -output lmder'];
%Get Intel Fortran Libraries (for LMDER build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lmdermex.c';
try
    eval([pre post])
    movefile(['lmder.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lmder','Error Compiling LM_DER / LM_DIF!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');