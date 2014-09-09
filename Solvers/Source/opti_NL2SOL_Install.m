%% NL2SOL Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile NL2SOL for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get NL2SOL
% NL2SOL is available in multiple variants, however the most recent is 
% included in the PORT library, available from:
% http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/port
% The library contains functions for a range of mathematical functions,
% however we will just be using the NL2SOL variants.

% 2) Compile NL2SOL + NL2SNO (DN2F, DN2G, DN2FB, DN2GB)
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

nl2path = 'C:\Solvers\port'; % FULL path to NL2SOL

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('NL2SOL',nl2path,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NL2SOL MEX file. Once you have completed all the
% above steps, simply run this file to compile NL2SOL! You MUST BE in the 
% base directory of OPTI!

clear nl2sol

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NL2SOL MEX FILE INSTALL\n\n');

%Get NL2SOL Libraries
post = [' -L' libdir ' -llibnl2sol -llibut -output nl2sol'];
%Get Intel Fortran Libraries (for NL2SOL build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nl2solmex.c';
try
    eval([pre post])
    movefile(['nl2sol.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nl2sol','Error Compiling NL2SOL!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');