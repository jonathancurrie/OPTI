%% MKL Trust Region NLS Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile the Intel Math Kernel Library (MKL) 
% Trust Region NLS function for use with MATLAB. NOTE you must NOT link the 
% threaded MKL libraries as the MATLAB callback function is not thread safe!

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get and Install Intel MKL
% http://software.intel.com/en-us/articles/intel-mkl/

% 6) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MKL TR NLS MEX file. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in the 
% base directory of OPTI!

clear mkltrnls

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL('seq'); %NOTE sequential only build!

fprintf('\n------------------------------------------------\n');
fprintf('MKL TR NLS MEX FILE INSTALL\n\n');

%Get MKL Libraries (for BLAS & Solver)
post = mkl_link;
%Common
post = [post ' -llibut -output mkltrnls'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims mkltrnls.c';
try
    eval([pre post])
    movefile(['mkltrnls.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:mkljac','Error Compiling MKL TR NLS!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
