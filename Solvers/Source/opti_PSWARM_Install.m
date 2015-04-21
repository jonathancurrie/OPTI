%% PSWARM Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile PSwarm for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get PSwarm
% PSwarm is available from http://www.norg.uminho.pt/aivaz/pswarm/. 
% Download PPSwarm_vxx.zip (C Version) and unzip to a suitable location.

% 2) Compile PSwarm
% The easiest way to compile PSwarm is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

pswmpath = 'C:\Solvers\PPSwarm_v1_5'; % FULL path to PSwarm 

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('PSwarm',pswmpath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the PSwarm MEX file. Once you have completed all the
% above steps, simply run this file to compile PSwarm! You MUST BE in the 
% base directory of OPTI!

clear pswarm

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('PSwarm MEX FILE INSTALL\n\n');

%Get PSwarm Libraries
post = [' -IInclude/Pswarm -L' libdir ' -llibpswarm -llibut -output pswarm'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims pswarmmex.c';
try
    eval([pre post])
    movefile(['pswarm.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:pswarm','Error Compiling PSwarm!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');