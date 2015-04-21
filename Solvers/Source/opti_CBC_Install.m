%% CBC Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Coin-Or Branch and Cut for use with 
% MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 0) Compile CLP as per the instructions in opti_CLP_Install.m

% 1) Get CBC
% CBC is available from https://projects.coin-or.org/Cbc. Download 
% the source.

% 2) Compile CBC
% The easiest way to compile CBC is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

cbcpath = 'C:\Solvers\Cbc-2.9.3\Cbc'; % FULL path to CBC

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('CBC',cbcpath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CBC MEX file. Once you have completed all 
% the above steps, simply run this file to compile CBC! You MUST BE in 
% the base directory of OPTI!

clear cbc

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('CBC MEX FILE INSTALL\n\n');

%Get CBC Libraries
post = [' -IInclude\Cbc -IInclude\Osi -IInclude\Cgl -L' libdir ' -llibcbc -llibcgl -llibut'];
%Get CLP and Osi libraries
post = [post ' -IInclude\Clp -IInclude\Coin -IInclude\Osi -IInclude\Cgl -IInclude\Cbc'];
post = [post ' -llibclp -llibcoinutils -llibosi -DCOIN_MSVS -output cbc'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims cbcmex.cpp';
try
    eval([pre post])
    movefile(['cbc.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:cbc','Error Compiling CBC!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
