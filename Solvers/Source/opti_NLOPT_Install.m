%% NLOPT Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile NonLinear OPTimization (NLOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 SP1 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get NLOPT
% NLOPT is available from http://ab-initio.mit.edu/wiki/index.php/NLopt. 
% Download the source.

% 2) Compile NLOPT
% The easiest way to compile NLOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):

nloptpath = 'C:\Solvers\nlopt-2.4.2'; % FULL path to NLOPT

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
%opti_VSBuild('NLOPT',nloptpath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NLOPT MEX file. Once you have completed all 
% the above steps, simply run this file to compile NLOPT! You MUST BE in 
% the base directory of OPTI!

clear nlopt

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NLOPT MEX FILE INSTALL\n\n');

%Get NLOPT Libraries
post = [' -IInclude\Nlopt -Inlopt\Include -L' libdir ' -llibnlopt -llibut -output nlopt'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nloptmex.c';
try
    eval([pre post])
    movefile(['nlopt.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nlopt','Error Compiling NLOPT!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
