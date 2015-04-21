%% LP_SOLVE Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile LP_SOLVE for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get LP_SOLVE
% LP_SOLVE is available from
% http://sourceforge.net/projects/lpsolve/files/lpsolve/. Download the
% source (lp_solve_5.5.2.0_source.tar.gz) or later version.

% 2) Compile LP_SOLVE
% The easiest way to compile LP_SOLVE is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:

lpspath = 'C:\Solvers\lp_solve_5.5'; % FULL path to LP_SOLVE

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
% opti_VSBuild('LPSOLVE',lpspath,cd,'VS2013');

% 3) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LP_SOLVE MEX file. Once you have completed all 
% the above steps, simply run this file to compile LP_SOLVE! You MUST BE in 
% the base directory of OPTI!

clear lp_solve

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LP_SOLVE MEX FILE INSTALL\n\n');

%Get LP_SOLVE Libraries
post = [' -IInclude\Lpsolve -Ilpsolve/Include -L' libdir ' -lliblpsolve -DMATLAB -DWIN32 -DLPSOLVEAPIFROMLIB -output lp_solve'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lpsolve/lpsolve.c lpsolve/matlab.c';
try
    eval([pre post])
    movefile(['lp_solve.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lpsolve','Error Compiling LP_SOLVE!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');