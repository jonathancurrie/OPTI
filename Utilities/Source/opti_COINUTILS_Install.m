%% CoinUtils Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile CoinUtils for use with MATLAB.

% My build platform:
% - Windows 7 SP1 x64
% - Visual Studio 2010

% 1) Compile Clp as per opti_CLP_Install.m in Solvers/Source

% 2) Compile the MEX Files
% The code below will automatically include all required libraries and
% directories to build the MEX files. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!

clear coinR coinW

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
solverdir = '..\..\Solvers\Source\';

fprintf('\n------------------------------------------------\n');
fprintf('COINUTILS MEX FILE INSTALL\n\n');

%Get COIN Libraries (NOTE I have called CoinUtils with GMPL libCoinUtilsGMPL)!!
post = [' -I' solverdir 'Include\Coin -I' solverdir 'Include\Glpk -L' solverdir libdir];
post = [post ' -llibcoinutilsgmpl -lglpk -DCOIN_MSVS '];
% post = [post ' -L' solverdir libdir ' -lglpk -I' solverdir '\Include\ '];

%CD to Source Directory
cdir = cd;
cd 'Utilities/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims coinR.cpp';
try
    eval([pre post])
     movefile(['coinR.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:coin','Error Compiling COINUTILS Read!\n%s',ME.message);
end
%Compile & Move
pre = 'mex -v -largeArrayDims coinW.cpp';
try
    eval([pre post])
     movefile(['coinW.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:coin','Error Compiling COINUTILS Write!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
