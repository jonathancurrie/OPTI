%% CLP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)
clc
% This file will help you compile Coin-Or Linear Programming for use with 
% MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get CLP
% CLP is available from http://www.coin-or.org/projects/Clp.xml. Download 
% the source. Unzip to a suitable temporary location.

% 2) Compile CLP & COIN Utils
% The easiest way to compile CLP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required paths on your computer:

clppath = 'C:\Solvers\Clp-1.16.5\Clp'; % FULL path to CLP
glpkpath = 'C:\Solvers\glpk-4.48'; % FULL path to GLPK(or leave blank [MAX VER 4.48])

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('CLP',{clppath,glpkpath},cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CLP MEX file. Once you have completed all 
% the above steps, simply run this file to compile CLP! You MUST BE in 
% the base directory of OPTI!

clear clp

%Enable for Aboca Build (must compile using Intel C++ & VS2012 Linker)
haveABC = false;

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('CLP MEX FILE INSTALL\n\n');

%Get CLP Libraries
post = [' -IInclude/Clp -IInclude/Coin -L' libdir ' -llibCoinUtils -llibut -llibosi -DCOIN_MSVS -output clp'];
%Get Optional Aboca
if(haveABC)
    post = [post ' -llibclpabc -DINTEL_COMPILER -DCLP_HAS_ABC=4 -D__BYTE_ORDER=__LITTLE_ENDIAN']; 
else
    post = [post ' -llibclp'];
end

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims clpmex.cpp';
try
    eval([pre post])
     movefile(['clp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:clp','Error Compiling CLP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
