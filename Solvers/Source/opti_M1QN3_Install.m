%% M1QN3 Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile M1QN3 for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get M1QN3
% M1QN3 is available from 
% https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html
% Download and Unzip the folder, we will create a Visual Studio project below.

% 2) Compile M1QN3
% The easiest way to compile M1QN3 is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

m1qpath = 'C:\Solvers\m1qn3-3.3-distrib'; % FULL path to M1QN3

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('M1QN3',m1qpath,cd);

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NL2SOL MEX file. Once you have completed all the
% above steps, simply run this file to compile M1QN3! You MUST BE in the 
% base directory of OPTI!

clear m1qn3

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('M1QN3 MEX FILE INSTALL\n\n');

%Get M1QN3 Libraries
post = [' -L' libdir ' -llibm1qn3 -llibut -output m1qn3'];
%Get Intel Fortran Libraries (for M1QN3 build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims m1qn3mex.c';
try
    eval([pre post])
    movefile(['m1qn3.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:m1qn3','Error Compiling M1QN3!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');