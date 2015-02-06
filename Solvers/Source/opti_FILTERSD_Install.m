%% FILTERSD Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile FILTER SD for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% 1) Get FILTER SD
% FilterSD is available from http://www.coin-or.org/download/source/filterSD/
% Download the .zip file for use on Windows.

% 2) Compile FILTERSD + FILTERSDSP
% The easiest way to compile FILTERSD is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

fsdpath = 'C:\Solvers\filterSD-1.0.0'; % FULL path to FilterSD

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('FilterSD',fsdpath,cd);

% 3) FILTERSD MEX Interface
% The FILTERSD MEX Interface is a simple MEX interface I wrote to use
% FILTERSD.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the FILTERSD MEX file. Once you have completed all the
% above steps, simply run this file to compile FILTERSD! You MUST BE in the 
% base directory of OPTI!

clear filtersd filtersdsp

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('FILTERSD MEX FILE INSTALL\n\n');

%Get FILTERSD Libraries
post = [' -L' libdir ' -llibfilterSD -llibut'];
%Get Intel Fortran Libraries (for ifort build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common Args
post = [post ' -output filtersd'];   

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims filtersdmex.c';
try
    eval([pre post])
    movefile(['filtersd.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:filtersd','Error Compiling FILTERSD!\n%s',ME.message);
end

%Compile & Move Sparse Version
post = regexprep(post,'-llibfilterSD','-llibfilterSDsp');
post = regexprep(post,'-output filtersd','-output filtersdsp');
try
    eval([pre post ' -DSPARSEVER'])
    movefile(['filtersdsp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:filtersd','Error Compiling FILTERSD SPARSE!\n%s',ME.message);
end

cd(cdir);
fprintf('------------------------------------------------\n');