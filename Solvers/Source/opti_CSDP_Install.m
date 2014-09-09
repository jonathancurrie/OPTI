%% CSDP Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile CSDP for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get CSDP
% CSDP is available from https://projects.coin-or.org/Csdp/. Download
% it and unzip the directory.

% 2) Compile CSDP
% The easiest way to compile CSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

csdppath = 'C:\Solvers\CSDP 6.2beta'; % FULL path to CSDP

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('CSDP',csdppath,cd);

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile CSDP! You MUST BE in the 
% base directory of OPTI!

clear csdp

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('CSDP MEX FILE INSTALL\n\n');

%Get Libraries
post = [' -IInclude/Csdp -L' libdir ' -llibcsdp -llibut -output csdp'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];


%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move (NOTE don't link against default VC++ OpenMP Lib as requires DLLs)
pre = 'mex -v -largeArrayDims LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:vcomp.lib" -DNOSHORTS csdpmex.c';
try
    eval([pre post])
    movefile(['csdp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:csdp','Error Compiling CSDP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');