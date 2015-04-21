%% IPOPT Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Interior Point OPTimizer (IPOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% NOTE - From OPTI v1.71 IPOPT is now dynamically linked against the
% MathWorks supplied libmwma57.dll (HSL MA57). This step is optional
% as we are also compiling MUMPS. Alternatively you can skip MUMPS, 
% and just use MA57! Also be aware libmwma57.dll does not play well
% on unconstrained problems due in part to missing MeTiS, thus the 
% ma57 pivot order option is overidden automatically.

% To recompile you will need to get / do the following:

% 1) Get IPOPT
% IPOPT is available from http://www.coin-or.org/download/source/Ipopt/.
% Download the latest version of the source.

% 2) Compile IPOPT
% The easiest way to compile IPOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

ipoptpath = 'C:\Solvers\Ipopt-3.12.1\Ipopt'; %FULL path to IPOPT
mumpspath = 'C:\Solvers\MUMPS_4.10.0'; %FULL path to MUMPS (or leave blank to skip linking MUMPS)

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
opti_VSBuild('IPOPT',{ipoptpath,mumpspath},cd,'VS2013');

% 3) IF LINKING MUMPS -> Complete MUMPS Compilation
% Complete opti_MUMPS_Install.m in Solvers/Source IF REQUIRED.

% 4) libmwma57 Compilation
% As MathWorks only supply a dll, I have manually created an import library
% for libmwma57.dll, so we can link against it. The project and source for
% this is included in libmwma57.zip. Simply open up the project, compile
% libmwma57.lib for your system (Win32/x64) and copy to 
%   Solvers/Source/lib/win32 or win64

% 5) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the IPOPT MEX file. Once you have completed all the
% above steps, simply run this file to compile IPOPT! You MUST BE in the 
% base directory of OPTI!

clear ipopt

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link,mkl_inc,mkl_lib,mkl_cmplr] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
% Compilation Options (THESE MUST MATCH THE OPTIONS IN opti_VSBuild.m)
haveMA57 = 1;
havePARDISO = 'MKL'; %MKL or BASEL version of PARDISO
haveLinearSolverLoader = 0; %Dynamically Load HSL's Linear Algebra Libraries
netlibBLAS = 0; %Link against deterministic (but slower) Netlib BLAS + LAPACK (requires you to compile them)

fprintf('\n------------------------------------------------\n');
fprintf('IPOPT MEX FILE INSTALL\n\n');

%Get IPOPT Libraries
post = ' -Iipopt/Include -IInclude/Ipopt -IInclude/BuildTools';
%Get MUMPS Libraries
if(~isempty(mumpspath))
    post = [post ' -IInclude\Mumps -llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
end
if(~haveLinearSolverLoader)
    %Optionally Add MA57
    if(haveMA57), post = [post ' -llibmwma57 -DhaveMA57']; end 
    %Optionally Add PARDISO
    switch(havePARDISO)
        case 'MKL'
            if(netlibBLAS), error('Cannot currently compile MKL_PARDISO with NETLIB BLAS'); end     
            post = [post ' -DhaveMKLPARDISO'];
        case 'BASEL'
            switch(computer)
                case 'PCWIN'
                    post = [post ' -llibpardiso412-WIN-X86 -DhavePARDISO'];
                case 'PCWIN64'
                    post = [post ' -llibpardiso412-WIN-X86-64 -DhavePARDISO'];
            end
    end
    post = [post ' -L' libdir ' -llibipopt'];
else
    post = [post ' -L' libdir ' -llibipoptdyn -llibma28part -DHAVE_LINEARSOLVERLOADER'];
end
%Select BLAS+LAPACK Libraries
if(netlibBLAS)
    post = [post ' -L"' mkl_cmplr '" -L"' cd '" -llapack -lblas ' mkl_for_link];
else
    %Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
    post = [post mkl_link mkl_for_link ' -DhaveMKL'];
end
%Common Args
post = [post ' -DIPOPT_BUILD -output ipopt'];   

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = ['mex -v -largeArrayDims ipopt/matlabexception.cpp ipopt/matlabfunctionhandle.cpp ipopt/matlabjournal.cpp '...
       'ipopt/iterate.cpp ipopt/ipoptoptions.cpp ipopt/options.cpp ipopt/sparsematrix.cpp ipopt/callbackfunctions.cpp '...
       'ipopt/matlabinfo.cpp ipopt/matlabprogram.cpp ipopt/ipopt.cpp'];
try
    eval([pre post])
    movefile(['ipopt.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ipopt','Error Compiling IPOPT!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');