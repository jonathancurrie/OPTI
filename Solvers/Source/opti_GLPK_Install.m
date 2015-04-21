%% GLPK Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile GNU Linear Programming Kit (GLPK) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get GLPK
% GLPK is available from http://www.gnu.org/software/glpk/glpk.html. We 
% will create VS projects below.

% 2) Compile GLPK
% The easiest way to compile GLPK is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

glpkpath = 'C:\Solvers\glpk-4.48'; % FULL path to GLPK [MAX VER 4.48]

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
%opti_VSBuild('GLPK',glpkpath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the GLPK MEX file. Once you have completed all the
% above steps, simply run this file to compile GLPK! You MUST BE in the 
% base directory of OPTI!

clear glpk

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('GLPK MEX FILE INSTALL\n\n');

%Get GLPK Libraries
post = [' -IInclude/Glpk -L' libdir ' -llibglpk -llibut -output glpk'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims glpkcc.cpp';
try
    eval([pre post])
    movefile(['glpk.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:glpk','Error Compiling GLPK!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');



