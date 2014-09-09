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

nmdpath = 'C:\Solvers\NOMAD-3.6.2'; % FULL path to NOMAD

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('Nomad',nmdpath,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NOMAD MEX file. Once you have completed all the
% above steps, simply run this file to compile NOMAD! You MUST BE in the 
% base directory of OPTI!

clear nomad

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NOMAD MEX FILE INSTALL\n\n');

%Get NOMAD Libraries
post = [' -IInclude/Nomad -L' libdir ' -llibnomad -llibut -output nomad -DOPTI_VERSION'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nomadmex.cpp';
try
    eval([pre post])
    movefile(['nomad.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nomad','Error Compiling NOMAD!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
