%% OPTI MEX Install Script

% This file enables me to quickly rebuild all solvers - it is not intended
% to be called by the user!

% My build platform:
% - Windows 10 x64
% - Visual Studio 2019
% - Intel Compiler XE (C++ & FORTRAN)
% - Intel Math Kernel Library (oneAPI)

% YOU MUST BE IN THE BASE DIRECTORY OF THE OPTI TOOLBOX!

addpath([cd '/Solvers/Source'])
addpath([cd '/Utilities/Source'])

% %% -- Intel C++ Solvers --
% clc
% fprintf(2,'Please Specify Intel C++ as your Compiler with VS2015 Linker (Option 3, Don''t Look for Installed Compilers)...\n\n');
% mex -setup
% %% CLP
% opti_CLP_Install


%% -- VS2019 Solvers --
clc
fprintf(2,'Please Specify Visual Studio 2019 as your Compiler...\n\n');
mex -setup c++
%% CLP
opti_CLP_Install

%% BONMIN
opti_BONMIN_Install

%% CBC
opti_CBC_Install

%% CSDP 
clc
% NOTE: Doesn't link under VS2019? Used 2015 and OK...
fprintf(2,'Please Specify Visual Studio 2015/2017 as your Compiler...\n\n');
mex -setup c++
%%
opti_CSDP_Install
%%
fprintf(2,'Please Specify Visual Studio 2019 as your Compiler...\n\n');
mex -setup c++

%% DSDP
opti_DSDP_Install

%% FILTERSD
opti_FILTERSD_Install

%% GLPK
opti_GLPK_Install

%% GSL
opti_GSL_Install

%% IPOPT
opti_IPOPT_Install

%% L-BFGS-B
opti_LBFGSB_Install

%% LEVMAR
opti_LEVMAR_Install

%% LP_SOLVE
opti_LPSOLVE_Install

%% M1QN3
opti_M1QN3_Install

%% MINPACK
opti_MINPACK_Install

%% MKL TR NLS
opti_MKLTRNLS_Install

%% MUMPS
opti_MUMPS_Install

%% NL2SOL
opti_NL2SOL_Install

%% NLOPT
opti_NLOPT_Install

%% NOMAD
opti_NOMAD_Install

%% OOQP
opti_OOQP_Install

%% PSWARM
opti_PSWARM_Install

%% SCIP
opti_SCIP_Install


%% -- Utilities --
%% MKL JAC
opti_MKLJAC_Install

%% CoinUtils
opti_COINUTILS_Install

%% AMPL
opti_ASL_Install

%% RMATHLIB
opti_RMathlib_Install


%% Remove Paths
rmpath([cd '/Solvers/Source'])
rmpath([cd '/Utilities/Source'])