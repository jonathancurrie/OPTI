
%% Visual Studio Builder Commands
path = 'C:\Users\Jonathan Currie\Desktop\Ipopt-3.11.3\Ipopt-3.11.3\Ipopt';%'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
mumpspath = 'C:\Users\Jonathan Currie\Documents\AUT\Code and Models\Optimization Solvers\MUMPS_4.10.0\MUMPS_4.10.0';%'full path to MUMPS here'; %e.g. 'C:\Solvers\MUMPS'
havePardiso = '';%'MKL'; %'BASEL' or 'MKL' %OPTIONALLY LINK PARDISO (CHOOSE VERSION TOO)
haveMA57 = 0;
sdir = [path '\src'];
hdrs = {[mumpspath '\include'], [mumpspath '\libseq'], [path '\..\BuildTools\headers'], [path '\..\ThirdParty\HSL']};
name = 'libipopt';
opts = [];
opts.exPP = {'IPOPT_BUILD','COIN_HAS_METIS','COIN_HAS_MUMPS','_CRT_SECURE_NO_WARNINGS','HAVE_LINEARSOLVERLOADER','HAVE_WINDOWS_H','SHAREDLIBEXT="dll"','COINHSL_HSL2013'};
switch(havePardiso)
    case 'BASEL'
        opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_PARALLEL'];
    case 'MKL'
        opts.exPP = [opts.exPP,'HAVE_PARDISO','HAVE_PARDISO_OLDINTERFACE','MKL_PARDISO'];
end
if(haveMA57)
    opts.exPP = [opts.exPP,'FUNNY_MA57_FINT'];
end
opts.exFolder = {'Apps','Algorithm\Inexact'}; %,'contrib\LinearSolverLoader'
opts.exclude = 'AmplTNLP.cpp';
opts.charset = 'multibyte';
VS_WriteProj(sdir,name,hdrs,opts)

%%
clear ipopt

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link,mkl_inc,mkl_lib,mkl_cmplr] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
mumpsdir = '..\..\mumps\Source\';
% Compilation Options
haveMA57 = 0;
havePARDISO = 0; %1 MKL, 2 BASEL
netlibBLAS = 0;

fprintf('\n------------------------------------------------\n');
fprintf('IPOPT MEX FILE INSTALL\n\n');

%Get IPOPT Libraries
post = ' -IInclude -IInclude/Common -IInclude/BuildTools -IInclude/Interfaces -IInclude/LinAlg -IInclude/LinAlg/TMatrices -IInclude/Algorithm ';
post = [post ' -L' libdir ' -llibIpopt'];
%GET HSL Libraries
post = [post ' -llibma28part'];
%Get MUMPS Libraries
post = [post ' -I' mumpsdir '\Include -L' mumpsdir libdir '-llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
%Optionally Add MA57
if(haveMA57), post = [post ' -llibmwma57 -DhaveMA57']; end
%Optionally Add PARDISO
switch(havePARDISO)
    case 1 %MKL
        if(netlibBLAS), error('Cannot currently compile MKL_PARDISO with NETLIB BLAS'); end     
        post = [post ' -DhaveMKLPARDISO'];
    case 2 %BASEL
        switch(computer)
            case 'PCWIN'
                post = [post ' -llibpardiso412-WIN-X86 -DhavePARDISO'];
            case 'PCWIN64'
                post = [post ' -llibpardiso412-WIN-X86-64 -DhavePARDISO'];
        end
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
cd 'Solvers/ipopt/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims matlabexception.cpp matlabfunctionhandle.cpp matlabjournal.cpp iterate.cpp ipoptoptions.cpp options.cpp sparsematrix.cpp callbackfunctions.cpp matlabinfo.cpp matlabprogram.cpp ipopt.cpp';
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


%% ALL NONLINEAR VERSION
clc
clear
% Objective & Gradient
 fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
 grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];
 
 % Constraints
A = [-1 1; 1 1];
rl = [-Inf;5];
ru = [-1;5];
lb = [0;0]; ub = [4;4];

% Nonlinear Constraint, Jacobian & Structure
 nlcon = @(x) A*x;
 nljac = @(x) sparse(A);
 jacstr = @() sparse(double(A~=0));        

% Starting Guess
 x0 = [2;2];

% Build Function Structure
 funcs.objective = fun;
 funcs.gradient = grad;
 funcs.constraints = nlcon;
 funcs.jacobian = nljac;
 funcs.jacobianstructure = jacstr;

% Build Options Structure
    opts = [];
 opts.lb = lb;
 opts.ub = ub;
 opts.cl = rl;
 opts.cu = ru;
 opts.ipopt.linear_solver = 'ma57';
 opts.ipopt.nlp_scaling_method = 'equilibration-based';
 opts.ipopt.hessian_approximation = 'limited-memory';
 opts.ipopt.dependency_detector = 'ma28';

% Call IPOPT
[x,output] = ipopt(x0,funcs,opts)