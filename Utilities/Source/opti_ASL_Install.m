%% ASL Install for OPTI Toolbox
% Supplied binaries are built from Netlib's AMPL Solver Library Interface

%   Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile AMPL Solver Library (ASL) for use with 
% MATLAB.

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% 1) Get AMPL Solver Library
% The generic NL reader for AMPL is available free from Netlib 
% (http://www.netlib.org/ampl/solvers/). You will need to download all .c
% and .h as well as .hd files. Note this is not the AMPL engine
% (www.ampl.com) which is a commerical product, but code to allow people to
% connect their solvers to AMPL. Alternatively send a blank email to 
% "netlib@netlib.org" with "send all from ampl/solvers" as the subject to 
% retrieve all files.

% 2) Compile ASL
% The easiest way to compile ASL is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

aslpath = 'C:\Solvers\ASL'; % FULL path to ASL

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('ASL',aslpath,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MEX file. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!

clear asl

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('AMPL MEX FILE INSTALL\n\n');

post = [' -IInclude/Asl -L' libdir ' -llibasl -DNO_STDIO1 -output asl'];

%CD to Source Directory
cdir = cd;
cd 'Utilities/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims amplmex.c';
try
    eval([pre post])
     movefile(['asl.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ampl','Error Compiling AMPL!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
