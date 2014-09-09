%% OOQP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Objective Orientated Quadratic
% Programming (OOQP) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get OOQP
% OOQP is available from http://pages.cs.wisc.edu/~swright/ooqp/. You will 
% need to register before you can download. We will create the VS project 
% below.

% 2) Compile OOQP
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

ooqppath = 'C:\Solvers\OOQP-0.99.22'; % FULL path to OOQP

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('OOQP',ooqppath,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the OOQP MEX file. Once you have completed all the
% above steps, simply run this file to compile OOQP! You MUST BE in the 
% base directory of OPTI!

clear ooqp

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('OOQP MEX FILE INSTALL\n\n');

%Get OOQP Libraries
post = [' -IInclude/Ooqp -L' libdir ' -llibooqp -llibut -output ooqp'];
%Get MA57 Library
post = [post ' -llibmwma57 -DHAVE_MA57'];
%Get Intel Fortran Libraries (for MA27 build if required) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link ' -DHAVE_PARDISO'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims ooqpmex.cpp';
try
    eval([pre post])
    movefile(['ooqp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ooqp','Error Compiling OOQP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
