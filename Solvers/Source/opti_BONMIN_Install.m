%% BONMIN Install for OPTI Toolbox
% Copyright (C) 2014 Jonathan Currie (I2C2)

% This file will help you compile Basic Open-source Nonlinear Mixed INteger
% programming (BONMIN) for use with MATLAB. 

% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 0) Complete Compilation as per OPTI instructions for CLP, CBC, MUMPS, and
% IPOPT, in that order.

% 1) Get BONMIN
% BONMIN is available from http://www.coin-or.org/Bonmin/. Download 
% the source. 

% 2) Compile BONMIN
% The easiest way to compile BONMIN is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

bminpath = 'C:\Solvers\Bonmin-1.8.1\Bonmin'; %FULL path to BONMIN
ipoptpath = 'C:\Solvers\Ipopt-3.12.1\Ipopt'; %FULL path to IPOPT
mumpspath = 'C:\Solvers\MUMPS_4.10.0'; %FULL path to MUMPS (or leave blank to skip linking MUMPS)
cbcpath = 'C:\Solvers\Cbc-2.9.3\Cbc'; % FULL path to CBC
clppath = 'C:\Solvers\Clp-1.16.5\Clp'; % FULL path to CLP

%Build VS Solution & Compile Solver Libraries (Win32 + Win64)
%opti_VSBuild('BONMIN',{bminpath,ipoptpath,mumpspath,cbcpath,clppath},cd,'VS2013');

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the BONMIN MEX file. Once you have completed all 
% the above steps, simply run this file to compile BONMIN! You MUST BE in 
% the base directory of OPTI!

clear bonmin 

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Compilation Options (THESE MUST MATCH THE OPTIONS IN opti_VSBuild.m)
haveMA57 = 1;
havePARDISO = 'MKL'; %MKL or BASEL

fprintf('\n------------------------------------------------\n');
fprintf('BONMIN MEX FILE INSTALL\n\n');

%Get Bonmin MEX Includes
post = ' -Ibonmin\Include -IInclude\BuildTools ';
%Get CBC & CGL Libraries
post = [post ' -IInclude\Cbc -IInclude\Osi -IInclude\Cgl -llibCbc -llibCgl -llibut'];
%Get CLP & COINUTILS libraries
post = [post ' -IInclude\Clp -IInclude\Coin -llibClp -llibOsi -llibCoinUtils'];
%Get MUMPS Libraries
if(~isempty(mumpspath))
    post = [post ' -IInclude\Mumps -llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
end
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Get IPOPT libraries
post = [post ' -IInclude\Ipopt -llibipopt'];
%Optionally Add MA57 and/or PARDISO
if(haveMA57), post = [post ' -llibmwma57 -DhaveMA57']; end 
%Optionally Add PARDISO
switch(havePARDISO)
    case 'MKL'  
        post = [post ' -DhaveMKLPARDISO'];
    case 'BASEL'
        switch(computer)
            case 'PCWIN'
                post = [post ' -llibpardiso412-WIN-X86 -DhavePARDISO'];
            case 'PCWIN64'
                post = [post ' -llibpardiso412-WIN-X86-64 -DhavePARDISO'];
        end
end
%Get BONMIN Includes
post = [post ' -IInclude/Bonmin ']; 
%Output Common
post = [post ' -L' libdir ' -DBONMIN_BUILD -DIPOPT_BUILD -llibbonmin -output bonmin'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Compile & Move
pre = ['mex -v -largeArrayDims bonmin/bonminmex.cpp bonmin/iterate.cpp bonmin/options.cpp bonmin/bonminoptions.cpp bonmin/matlabinfo.cpp '...
       'bonmin/callbackfunctions.cpp bonmin/matlabexception.cpp bonmin/matlabfunctionhandle.cpp bonmin/matlabprogram.cpp bonmin/matlabjournal.cpp bonmin/sparsematrix.cpp'];
try
    eval([pre post])
    movefile(['bonmin.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:bonmin','Error Compiling BONMIN!\n%s',ME.message);
end

% %BUILD CPLEX BONMIN VERSION
% if(haveCPLEX)
%     % Modify below function it it cannot find IBM ILOG CPLEX on your system
%     [CPLX_link] = opti_FindCplex();
%     post = [post CPLX_link ' -llibosicpx -DHAVE_CPLEX' ];
%     % Include BONMIN build with CPLEX
%     post = regexprep(post,' -llibbonmin -output bonmin',' -llibbonminCplex -output bonminCplex');
%     
%     %Compile & Move
%     pre = ['mex -v -largeArrayDims bonminmex.cpp iterate.cpp options.cpp bonminoptions.cpp matlabinfo.cpp '...
%            'callbackfunctions.cpp matlabexception.cpp matlabfunctionhandle.cpp matlabprogram.cpp matlabjournal.cpp sparsematrix.cpp'];
%     try
%         eval([pre post])
%         movefile(['bonminCplex.' mexext],'../','f')
%         fprintf('Done!\n');
%     catch ME
%         cd(cdir);
%         error('opti:bonmin','Error Compiling BONMIN CPLEX Version!\n%s',ME.message);
%     end
% end

cd(cdir);
fprintf('------------------------------------------------\n');
