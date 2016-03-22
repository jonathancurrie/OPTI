%% NOMAD Install for OPTI Toolbox

% This file will help you compile NOMAD for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from http://www.gerad.ca/NOMAD/PHP_Forms/Download.php.
% Complete the download form then download the latest version. Once you
% have installed NOMAD, locate the /src/ directory.

% 2) Compile NOMAD
% The easiest way to compile NOMAD is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% path = 'C:\Solvers\NOMAD-3.7.2'; % FULL path to NOMAD
% opti_VSBuild('Nomad',path);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NOMAD MEX file. Once you have completed all the
% above steps, simply run this file to compile NOMAD! You MUST BE in the 
% base directory of OPTI!

%MEX Interface Source Files
src = 'nomadmex.cpp';
%Include Directories
inc = 'Include/Nomad';
%Lib Names [static libraries to link against]
libs = 'libnomad';
%Options
opts = [];
opts.verb = false;
opts.pp = 'OPTI_VERSION';

%Compile
opti_solverMex('nomad',src,inc,libs,opts);