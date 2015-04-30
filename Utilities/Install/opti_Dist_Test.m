%% Pre-Distribution-Release Test File
% Variety of tests to check opti OK for distribution
clc
clear all

%Run Install Test
opti_Install_Test;

%Check MEX Interfaces all working
vers = checkSolver('ver');

%Run Benchmarks
optiBench('LP');
optiBench('BILP');
optiBench('MILP');
optiBench('QP');
optiBench('SDP');
optiBench('NLS');
optiBench('NLP');

%Run Demos
Basic_Usage;
MATLAB_Overloads;
LinearProgramming;
MixedInteger_LinearProgramming;
QuadraticProgramming;
Global_NonlinearProgramming;
NonlinearProgramming;
Differentiation_Examples;
FileIO_Examples;
AMPL_Examples;

%Run MATLAB Wrapper Soak Tests
test_probs_mw;
test_probs_mw;
test_probs_mw;

%Run MPS Test
test_mps;

%Run SOS Test
test_sos;

%QC Rest
test_new_qc;

%Run Row Bounds Test
test_rowcon

%Run NumDiff Test
test_numdiff;

%Run Global Test Set
test_global;

%Run iterfun Test
test_iterfunc;

%Run hessian Test
test_hessian;

%Run AMPL Soak Tests
test_ampl;
test_ampl;
test_ampl;

%Run Solver Tests
test_clp;
test_cbc;
test_filterSD;
test_nlopt;
test_ooqp;
test_scip_formal;

%Run linear constraint tests
test_ipopt_linearcon;
test_bonmin_linearcon;

%Run sdp tests
test_csdp;
test_dsdp;
test_sedumi;
test_sdpio;

%Run other tests
test_objc;
test_sym_diff;
test_sparse_snle;
test_rmathlib_formal;

%Plot tests
test_1dplots;
test_3dplots;

%Multisolve Tests
test_multisolve;
auto2fit_problems;

%Confidence Tests
test_conf_formal;
test_conf;
test_optifit;
test_nist_fit;

%Dynamic Optimization Tests
test_dyn_nls;
dnls_examples;

%Run General Test Set
test_probs;

clc
fprintf('\nAll Appears OK!\n\n');
