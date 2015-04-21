%% DSDP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile DSDP for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get DSDP
% DSDP is available from http://www.mcs.anl.gov/hs/software/DSDP/. Download
% it and unzip the directory.

% 2) Compile DSDP
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

dsdppath = 'C:\Solvers\DSDP5.8'; % FULL path to DSDP

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
%opti_VSBuild('DSDP',dsdppath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the DSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile DSDP! You MUST BE in the 
% base directory of OPTI!

clear dsdp

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('DSDP MEX FILE INSTALL\n\n');

%Get Libraries
post = [' -IInclude/DSDP -L' libdir ' -llibdsdp -llibut -output dsdp'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims dsdpmex.c';
try
    eval([pre post])
    movefile(['dsdp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:dsdp','Error Compiling DSDP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');