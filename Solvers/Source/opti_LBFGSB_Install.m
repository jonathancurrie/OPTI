%% L-BFGS-B Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Limited Memory Broyden-Fletcher-Goldfarb-
% Shanno Bounded Optimization (L-BFGS-B) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get L-BFGS-B
% L-BFGS-B is available from 
% http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html. We will create 
% VS projects below.

% 2) Compile L-BFGS-B
% The easiest way to compile L-BFGS-B is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

lbfgsbpath = 'C:\Solvers\Lbfgsb.3.0'; % FULL path to L-BFGS-B

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
opti_VSBuild('LBFGSB',lbfgsbpath,cd,'VS2013');
%%
% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the L-BFGS-B MEX file. Once you have completed all the
% above steps, simply run this file to compile L-BFGS-B! You MUST BE in the 
% base directory of OPTI!

clear lbfgsb

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('L-BFGS-B MEX FILE INSTALL\n\n');

%Get LBFGSB Libraries
post = [' -Ilbfgsb\Include -L' libdir ' -llibLBFGSB -llibut'];
%Get Intel Fortran Libraries (for LBFGSB build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lbfgsb/lbfgsb.cpp lbfgsb/lbfgsb_program.cpp lbfgsb/program.cpp';
try
    eval([pre post])
    movefile(['lbfgsb.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lbfgsb','Error Compiling L-BFGS-B!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');