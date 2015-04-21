%% LEVMAR Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Levenberg-Marquardt in C/C++ (LEVMAR) for 
% use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get LEVMAR
% LEVMAR is available from http://www.ics.forth.gr/~lourakis/levmar/. We 
% will create the VS project below.

% 2) Compile LEVMAR
% The easiest way to compile LEVMAR is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

levpath = 'C:\Solvers\levmar-2.6'; %FULL path to LEVMAR

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('LEVMAR',levpath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LEVMAR MEX file. Once you have completed all the
% above steps, simply run this file to compile LEVMAR! You MUST BE in the 
% base directory of OPTI!

clear levmar

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LEVMAR MEX FILE INSTALL\n\n');

%Get LEVMAR Libraries
post = [' -IInclude/Levmar -L' libdir ' -lliblevmar'];
%Get MKL Libraries (for BLAS & LAPACK)
post = [post mkl_link];
%Common outputs
post = [post ' -output levmar'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims levmarmex.c';
try
    eval([pre post])
    movefile(['levmar.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:levmar','Error Compiling LEVMAR!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
