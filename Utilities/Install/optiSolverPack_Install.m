%OPTI Solver Package Installation
clc
clear
% This file will help you compile two DLLs:

%optiSolverPackA containing:
% - CoinUtils
% - Mumps
% - Metis
% - Clp
% - Glpk
% - Cbc
% - IPOPT
% - Bonmin
% - OOQP
% - Nomad (FIX!!)

%optiSolverPackB containing:
% - LBFGSB
% - Levmar
% - M1QN3
% - HYBRJ
% - LMDER
% - NL2SOL
% - NLOPT
% - PSwarm (FIX!!)
% - Lp_solve
% - DSDP

% - RMathlib
% - CSDP
% - ASL

% In addition, it will create static libraries suitable for linking the
% following solvers
% - FilterSD (dense and sparse version)


% My build platform:
% - Windows 7 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get the following packages (always the source):
% BONMIN  - Available from http://www.coin-or.org/Bonmin/.
% GLPK - 
% Mumps - 
% Metis -

% LMDER+HYBRJ - Available from http://www.netlib.org/minpack/. (Download and combine lmder, lmdir, hybrd + hybrj + all dependencies into single directory)

% 2) Use Visual Studio Builder to Create the Project
% The easiest way to compile is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

%% Visual Studio Builder Commands

path = 'C:\Solvers';
bminpath = 'C:\Solvers\Bonmin-1.7.4'; %full path to BONMIN Package
glpkpath = 'C:\Solvers\glpk-4.48'; %full path to GLPK
mumpspath = 'C:\Solvers\MUMPS_4.10.0'; %full path to MUMPS
metispath = 'C:\Solvers\metis-4.0.3'; %full path to METIS
filsdpath = 'C:\Solvers\filterSD-1.0.0'; %full path to filterSD
lbfgspath = 'C:\Solvers\Lbfgsb.3.0'; %full path to L-BFGS-B
levmpath = 'C:\Solvers\levmar-2.6'; %full path to LEVMAR
m1qpath = 'C:\Solvers\m1qn3-3.3-distrib'; %full path to M1QN3
mnpkpath = 'C:\Solvers\minpack'; %full path to the MINPACK library (for HYBRJ + LMDER)
portpath = 'C:\Solvers\port'; %full path to the PORT library (for NL2SOL)
ooqppath = 'C:\Solvers\OOQP-0.99.22'; %full path to OOQP
pswmpath = 'C:\Solvers\PPSwarm_v1_5'; %full path to PSwarm
nloptpath = 'C:\Solvers\nlopt-2.4.2'; %full path to Nlopt
dsdppath = 'C:\Solvers\DSDP5.8'; %full path to DSDP
lpspath = 'C:\Solvers\lp_solve_5.5'; %full path to LP_SOLVE
nmdpath = 'C:\Solvers\NOMAD-3.6.2'; %full path to NOMAD
csdppath = ''; %full path to CSDP
aslpath = ''; %full path to ASL
rmathpath = ''; %full path to RMathlib

libsA = {'libbonmin.lib','libipopt.lib','libcbc.lib','libclp.lib','libcgl.lib','libcoinutils.lib','libosi.lib',...
           'libdmumps_c.lib','libdmumps_f.lib','libzmumps_c.lib','libzmumps_f.lib','libseq_c.lib','libseq_f.lib',...
           'libmetis.lib','libpord.lib','libglpk.lib','libooqp.lib','libnomad.lib'};
       
libsB = {'liblbfgsb.lib','liblevmar.lib','libm1qn3.lib','libminpack.lib','libnl2sol.lib','libpswarm.lib',...
         'libnlopt.lib','libdsdp.lib','liblpsolve.lib'};       

n = 1;
%OPTISOLVERPACKA
sdir = path;
hdrs = [];
pkgnameA = 'optiSolverPackA';
opts = [];
opts.empty = true;
opts.dll = true;
opts.def = [pkgnameA '.def'];
opts.linklib = [libsA 'libmwma57.lib'];
opts.linkpath = {'$(SolutionDir)'};
opts.mkllink = true;
opts.archname = true;
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=pkgnameA; VSPRJ(n).opts=opts; n = n + 1;
%OPTISOLVERPACKB
sdir = path;
hdrs = [];
pkgnameB = 'optiSolverPackB';
opts = [];
opts.empty = true;
opts.dll = true;
opts.def = [pkgnameB '.def'];
opts.linklib = libsB;
opts.linkpath = {'$(SolutionDir)'};
opts.mkllink = true;
opts.archname = true;
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=pkgnameB; VSPRJ(n).opts=opts; n = n + 1;
%BONMIN
%Comment asl define
str = fileread([bminpath '/Bonmin/src/Interfaces/config_default.h']);
if(isempty(strfind(str,'//#define COIN_HAS_ASL 1')))
    str = regexprep(str,'#define COIN_HAS_ASL 1','//#define COIN_HAS_ASL 1');
    fp = fopen([bminpath '/Bonmin/src/Interfaces/config_default.h'],'w');
    fprintf(fp,'%s',str); fclose(fp);
end
sdir = [bminpath '\Bonmin\src'];
hdrs = {[bminpath '\Cgl\src'],[bminpath '\Cbc\src'],[bminpath '\Clp\src'],[bminpath '\CoinUtils\src'],[bminpath '\Ipopt\src'],[bminpath '\Osi\src'], [bminpath '\BuildTools\headers']};
name = 'libbonmin';
opts = [];
opts.exPP = {'BONMIN_BUILD','IPOPT_BUILD','FUNNY_MA57_FINT','COINHSL_HAS_MA57','_CRT_SECURE_NO_WARNINGS',...
             'HAVE_PARDISO','HAVE_PARDISO_OLDINTERFACE','HAVE_PARDISO_MKL'};
opts.exclude = {'BonCurvatureEstimator.cpp','BonCurvBranchingSolver.cpp'};
opts.exFolder = {'Apps','Interfaces\Ampl','Interfaces\Filter','Algorithms\Ampl'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%CLP
clppath = [bminpath '\Clp'];
sdir = [clppath '\src'];
hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], [clppath '\..\BuildTools\headers']};
name = 'libclp';
opts = [];
opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
                'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
opts.exFilter = {'Abc*','CoinAbc*'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% CoinUtils
sdir = [clppath '\..\CoinUtils\src'];
hdrs = {[glpkpath '\src'], [clppath '\..\BuildTools\headers']};
name = 'libcoinutils';
opts = [];
opts.exPP = {'COIN_HAS_GLPK','COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% Osi
sdir = [clppath '\..\Osi\src'];
hdrs = {[clppath '\..\CoinUtils\src'] [clppath '\..\BuildTools\headers']};
name = 'libosi';
opts = [];
opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% CBC
cbcpath = [bminpath '\Cbc'];
sdir = [cbcpath '\src'];
hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Cgl\src'], [cbcpath '\..\Clp\src'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
name = 'libcbc';
opts = [];
opts.exPP = {'CBC_BUILD','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
opts.exclude = {'unitTest.cpp','CoinSolve.cpp','CbcGeneric.cpp','CbcBranchBase.cpp'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% CGL
sdir = [cbcpath '\..\Cgl\src'];
hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Clp\src\'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
name = 'libcgl';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};          
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% IPOPT 
%Comment hsl & asl defines
str = fileread([bminpath '/Ipopt/src/Common/config_default.h']);
if(isempty(strfind(str,'//#define COIN_HAS_ASL 1')))
    str = regexprep(str,'#define COIN_HAS_ASL 1','//#define COIN_HAS_ASL 1');
    str = regexprep(str,'#define COIN_HAS_HSL 1','//#define COIN_HAS_HSL 1');
    fp = fopen([bminpath '/Ipopt/src/Common/config_default.h'],'w');
    fprintf(fp,'%s',str); fclose(fp);
end
ippath = [bminpath '\Ipopt'];
sdir = [ippath '\src'];
hdrs = {[mumpspath '\include'], [mumpspath '\libseq'], [ippath '\..\BuildTools\headers']};
name = 'libipopt';
opts = [];
opts.exPP = {'IPOPT_BUILD','FUNNY_MA57_FINT','COINHSL_HAS_MA57','COIN_HAS_METIS','COIN_HAS_MUMPS',...
             '_CRT_SECURE_NO_WARNINGS','HAVE_PARDISO','HAVE_PARDISO_OLDINTERFACE','HAVE_PARDISO_MKL'};
opts.exFolder = {'Apps','Algorithm\Inexact','contrib\LinearSolverLoader'};
opts.exclude = 'AmplTNLP.cpp';
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;    
% DMUMPS_C
sdir = [mumpspath '\src'];
hdrs = {[mumpspath '\include'],[mumpspath '\PORD\include']};
name = 'libdmumps_c';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_d','pord','metis'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% DMUMPS_F
sdir = [mumpspath '\src'];
hdrs = {[mumpspath '\include'],[mumpspath '\libseq']};
name = 'libdmumps_f';
opts = [];
opts.cpp = false;
opts.exPP = {'pord','metis'};
opts.exFilter = {'cmumps*','smumps*','zmumps*'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% ZMUMPS_C 
sdir = [mumpspath '\src'];
hdrs = {[mumpspath '\include'],[mumpspath '\PORD\include']};
name = 'libzmumps_c';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_z','pord','metis'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% ZMUMPS_F 
sdir = [mumpspath '\src'];
hdrs = {[mumpspath '\include'],[mumpspath '\libseq']};
name = 'libzmumps_f';
opts = [];
opts.cpp = false;
opts.exPP = {'pord','metis'};
opts.exFilter = {'cmumps*','smumps*','dmumps*'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% LIBSEQ_C
sdir = [mumpspath '\libseq'];
name = 'libseq_c';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% LIBSEQ_F
sdir = [mumpspath '\libseq'];
name = 'libseq_f';
opts = [];       
opts.cpp = false;
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% PORD
sdir = [mumpspath '\PORD\lib'];
hdrs = [mumpspath '\PORD\include'];
name = 'libpord';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% METIS
%Comment strings
str = fileread([metispath '/Lib/metis.h']);
str = regexprep(str,'#include <strings.h>','//#include <strings.h>');
fp = fopen([metispath '/Lib/metis.h'],'w');
fprintf(fp,'%s',str);
fclose(fp);
sdir = [metispath '\Lib'];
hdrs = [];
name = 'libmetis';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','__STDC__','__VC__'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% GLPK
sdir = [glpkpath '\src'];
hdrs = [];
name = 'libglpk';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% FILTERSD (dense)
sdir = filsdpath;
name = 'libfilterSD';
opts = [];       
opts.cpp = false;
opts.include = {'checkd.f','filterSD.f','glcpd.f','l1sold.f','denseL.f','denseA.f','util.f'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% FILTERSD (sparse)
sdir = filsdpath;
name = 'libfilterSDsp';
opts = [];       
opts.cpp = false;
opts.include = {'checkd.f','filterSD.f','glcpd.f','l1sold.f','schurQR.f','sparseA.f','util.f'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% L-BFGS-B
sdir = lbfgspath;
name = 'liblbfgsb';
opts = [];       
opts.cpp = false;
opts.include = {'lbfgsb.f','timer.f','linpack.f'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% LEVMAR
sdir = levmpath;
name = 'liblevmar';
opts = [];       
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
opts.exclude = {'expfit.c','Axb_core.c','lm_core.c','lmbc_core.c','lmblec_core.c','lmbleic_core.c','lmlec_core.c','misc_core.c'};
opts.exFolder = {'matlab'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% M1QN3
sdir = [m1qpath '\src'];
name = 'libm1qn3';
opts = [];       
opts.cpp = false;
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% MINPACK (LMDER + HYBRJ)
sdir = mnpkpath;
name = 'libminpack';
opts = [];       
opts.cpp = false;
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% NL2SOL
%Fix machine constants
str = fileread([portpath '/Mach/d1mach.f']);
str = strrep(str,sprintf('\n      DATA SMALL(1),SMALL(2) /    1048576,          0 /'),sprintf('\nC      DATA SMALL(1),SMALL(2) /    1048576,          0 /'));
str = strrep(str,sprintf('\n      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /'),sprintf('\nC      DATA LARGE(1),LARGE(2) / 2146435071,         -1 /'));
str = strrep(str,sprintf('\n      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /'),sprintf('\nC      DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /'));
str = strrep(str,sprintf('\n      DATA DIVER(1),DIVER(2) / 1018167296,          0 /'),sprintf('\nC      DATA DIVER(1),DIVER(2) / 1018167296,          0 /'));
str = strrep(str,sprintf('\n      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/'),sprintf('\nC      DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/'));
str = strrep(str,'C      DATA SMALL(1),SMALL(2) /          0,    1048576 /','      DATA SMALL(1),SMALL(2) /          0,    1048576 /');
str = strrep(str,'C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /','      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /');
str = strrep(str,'C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /','      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /');
str = strrep(str,'C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /','      DATA DIVER(1),DIVER(2) /          0, 1018167296 /');
str = strrep(str,'C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/','      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/');
fp = fopen([portpath '/Mach/d1mach.f'],'w');
fprintf(fp,'%s',str); fclose(fp);
sdir = portpath;
name = 'libnl2sol';
opts = [];       
opts.cpp = false;
opts.include = {'dn2f.f','dn2g.f','dn2fb.f','dn2gb.f','dn2rdp.f','dv7scp.f','divset.f','Mach\d1mach.f','Mach\r1mach.f',...
                'Mach\i1mach.f','drn2g.f','drn2gb.f','dd7tpr.f','dd7upd.f','dg7itb.f','ditsum.f','dl7vml.f',...
                'dq7apl.f','dq7rad.f','dr7tvm.f','dv7cpy.f','dv2nrm.f','dc7vfn.f','dg7lit.f','dn2cvp.f',...
                'dn2lrd.f','da7sst.f','df7dhb.f','df7hes.f','dg7qsb.f','dg7qts.f','dl7itv.f','dl7ivm.f',...
                'dl7msb.f','dl7mst.f','dl7nvr.f','dl7sqr.f','dl7srt.f','dl7svn.f','dl7svx.f','dl7tsq.f',...
                'dl7tvm.f','do7prd.f','dparck.f','dq7rsh.f','dr7mdc.f','drldst.f','ds7dmp.f','ds7ipr.f',...
                'ds7lup.f','ds7lvm.f','dv2axy.f','dv7dfl.f','dv7ipr.f','dv7scl.f','dv7vmp.f','i7copy.f',...
                'i7mdcn.f','i7pnvr.f','i7shft.f','stopx.f','dh2rfa.f','dh2rfg.f','ds7bqn.f','fdump.f',...
                'seterr.f','dd7mlp.f','dv7shf.f','e9rint.f','eprint.f','i8save.f','s88fmt.f','sdump.f',...
                'stkdmp.f','a9rntc.f','a9rntd.f','a9rnti.f','a9rntl.f','a9rntr.f','frmatd.f','frmati.f',...
                'frmatr.f','i0tk00.f','u9dmp.f','i10wid.f','iceil.f','iflr.f'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%OOQP
%Fix couple problems
str = fileread([ooqppath '/src/Readers/hash.C']);
str = regexprep(str,'hash(','hashL(');
fp = fopen([ooqppath '/src/Readers/hash.C'],'w');
fprintf(fp,'%s',str); fclose(fp);
str = fileread([ooqppath '/src/QpSolvers/MehrotraSolver.C']);
if(isempty(strfind(str,'extern double gmu;')))
    str = regexprep(str,'double gmu;','extern double gmu;');
    fp = fopen([ooqppath '/src/QpSolvers/MehrotraSolver.C'],'w');
    fprintf(fp,'%s',str); fclose(fp);
end
sdir = [ooqppath '\src'];
name = 'libooqp';
opts = [];  
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
opts.mkl = true;
opts.exclude = {'QpGenSparseOblio.C','QpGenSparseSuperLu.C','OoqpPetscMonitor.C',...
                'QpGenSparseGondzioDriver.C','QpGenSparseMa57GondzioDriver.C','QpGenSparseMehrotraDriver.C',...
                'QpBoundDenseGondzioDriver.C','QpGenDenseGondzioDriver.C','SvmGondzioDriver.C','HuberGondzioDriver.C'};
opts.exFilter = {'QpBoundPetsc*'};
opts.exFolder = {'Ampl','CInterface','PetscLinearAlgebra','QpExample','Mex'};
opts.compileAsCpp = true;
opts.ExcepwCExtern = true;
%Copy over OPTI modifications
unzip([cd '/Solvers/Source/ooqp/Ma57Solver.zip'],[cd '/Solvers/Source/ooqp/Ma57']);
unzip([cd '/Solvers/Source/ooqp/PardisoSolver.zip'],[cd '/Solvers/Source/ooqp/Pardiso']);
pause(0.1); rehash;
copyfile([cd '/Solvers/Source/ooqp/Ma57/Ma57Solver.c'],[ooqppath '\src\LinearSolvers\Ma57Solver\Ma57Solver.c'],'f');
copyfile([cd '/Solvers/Source/ooqp/Ma57/Ma57Solver.h'],[ooqppath '\src\LinearSolvers\Ma57Solver\Ma57Solver.h'],'f');
if(~exist([ooqppath '\src\LinearSolvers\Pardiso'])), mkdir([ooqppath '\src\LinearSolvers\'],'Pardiso'); end
copyfile([cd '/Solvers/Source/ooqp/Pardiso/PardisoSolver.cpp'],[ooqppath '\src\LinearSolvers\Pardiso\PardisoSolver.cpp'],'f');
copyfile([cd '/Solvers/Source/ooqp/Pardiso/PardisoSolver.h'],[ooqppath '\src\LinearSolvers\Pardiso\PardisoSolver.h'],'f');
copyfile([cd '/Solvers/Source/ooqp/Pardiso/QpGenSparsePardiso.cpp'],[ooqppath '\src\QpGen\QpGenSparsePardiso.cpp'],'f');
copyfile([cd '/Solvers/Source/ooqp/Pardiso/QpGenSparsePardiso.h'],[ooqppath '\src\QpGen\QpGenSparsePardiso.h'],'f');
rmdir([cd '/Solvers/Source/ooqp/Pardiso'],'s');
rmdir([cd '/Solvers/Source/ooqp/Ma57'],'s');
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%PSWARM
sdir = pswmpath;
name = 'libpswarm';
opts = [];
opts.exPP = {'LINEAR','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
opts.exclude = {'showcache.c','user.c','pswarm_main.c','pswarm_py.c','pswarm_r.c','cache.c'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%NLOPT
%Copy config.h
copyfile([cd '/Solvers/Source/nlopt/config.h'],[nloptpath '\api\config.h'],'f');
sdir = nloptpath;
name = 'libnlopt';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
opts.exclude = {'testfuncs.c','tst.cc','tstc.c','testros.cc','prog.cc','redblack_test.c','DIRparallel.c'};
opts.exFolder = {'octave','swig','test'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%DSDP
sdir = [dsdppath '\src'];
hdrs = [dsdppath '\include'];
name = 'libdsdp';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%LPSOLVE
%Copy lusol.c + .h
copyfile([lpspath '\bfp\bfp_LUSOL\LUSOL\lusol.c'],[lpspath '\lusol.c'],'f');
copyfile([lpspath '\bfp\bfp_LUSOL\LUSOL\lusol.h'],[lpspath '\lusol.h'],'f');
%Copy modified myblas (to avoid hassles with mkl vs local blas)
unzip([cd '/Solvers/Source/lpsolve/myblas.zip'],[cd '/Solvers/Source/lpsolve/myblas']);
pause(0.1); rehash;
copyfile([cd '/Solvers/Source/lpsolve/myblas/myblas.c'],[lpspath '\shared\myblas.c'],'f');
rmdir([cd '/Solvers/Source/lpsolve/myblas'],'s');
sdir = lpspath;
hdrs = [lpspath '\bfp\bfp_LUSOL\LUSOL'];
name = 'liblpsolve';
opts = [];
opts.exPP = {'PARSER_LP','INVERSE_ACTIVE=INVERSE_LUSOL','RoleIsExternalInvEngine',...
             '_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','LoadableBlasLib=0'};
opts.exclude = {'lp_solveDLL.c'};
opts.exFilter = {'lp_BFP*'};
opts.exFolder = {'lpsolve55','lp_solve','demo','bfp\bfp_LUSOL\LUSOL'}; 
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%Nomad
%Increase max dim
str = fileread([nmdpath '/defines.hpp']);
if(isempty(strfind(str,'const int MAX_DIMENSION = 5000;')))
    str = regexprep(str,'const int MAX_DIMENSION = 1000;','const int MAX_DIMENSION = 5000;');
    fp = fopen([nmdpath '/defines.hpp'],'w');
    fprintf(fp,'%s',str); fclose(fp);
end
%Fix static method prob + function defined in header issue
str = fileread([nmdpath '/Eval_Point.hpp']); 
rep = 'void reset_tags_and_bbes ( void ) {_current_tag = 0;_current_bbe = 0;_current_sgte_bbe = 0;}';
rep1 = 'bool is_eval_ok ( void ) const { return (_eval_status == NOMAD::EVAL_OK); }';
rep2 = 'bool is_feasible ( const NOMAD::Double & h_min ) const';
if(~isempty(strfind(str,rep)))
    str = strrep(str,rep,'void reset_tags_and_bbes ( void );');
    str = strrep(str,rep1,'bool is_eval_ok ( void ) const;');
    str = strrep(str,rep2,'bool is_feasible ( const NOMAD::Double & h_min ) const;');
    str = strrep(str,'return ( _h.is_defined() && _h <= h_min );','');
    fp = fopen([nmdpath '/Eval_Point.hpp'],'w');
    fprintf(fp,'%s',str); fclose(fp);
    str = fileread([nmdpath '/Eval_Point.cpp']); 
    str = sprintf('%s\n\n%s',str,'void NOMAD::Eval_Point::reset_tags_and_bbes ( void ) {_current_tag = 0;_current_bbe = 0;_current_sgte_bbe = 0;}');
    str = sprintf('%s\n\n%s',str,'bool NOMAD::Eval_Point::is_eval_ok ( void ) const { return (_eval_status == NOMAD::EVAL_OK); }');
    str = sprintf('%s\n\n%s',str,'bool NOMAD::Eval_Point::is_feasible ( const NOMAD::Double & h_min ) const { return ( _h.is_defined() && _h <= h_min ); }');
    fp = fopen([nmdpath '/Eval_Point.cpp'],'w');
    fprintf(fp,'%s',str); fclose(fp);
end
sdir = nmdpath;
name = 'libnomad';
opts = [];
opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
opts.exclude = 'nomad.cpp';
VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%Write the Solution File
% VS_WriteSol(VSPRJ)


%% Open the Visual Studio Solution (optiSolverPackA)
% Build -> Configuration Manager -> Remove optiSolverPackA+B from Build for
% Release Builds

%Make the REQUIRED Changes below:

%NOMAD
%   a) In Eval_Point.hpp line 384 delete the empty curly braces

%FILTERSD
%   a) Replace filterSD.f and glcpd.f with the files located in the
%   following folder (only compatible with OPTI version):
%       OPTI/Solvers/Source/filterSD/filterSD_JCEdit.zip

%PSWARM (FIX ME)
%   a) A couple code changes are required in order to interface with OPTI:
%       - In pswarm.h add "int solveriters;" to the Stats structure as the last line.
%       - In pswarm.h change the define on line 62 to "#define SYS_RANDOM 1"
%       - In pswarm.c on line 666 add "stats.solveriters = iter;"
%       - Change the structure definitions in pswarm.c at line 66 to extern:
%           "extern struct Stats stats;"
%       - Add the following line under the above:
%           "extern struct Options opt;"
%       - Comment the default options structure in pswarm.c on line 70
%       (comment all lines for this structure).



%Right click the solution, build.

%% Generate Module Definition File (MATLAB Arch Dependent)
clc
cdir = cd;
vcdir = 'C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC';
if(strcmp(computer,'PCWIN'))
%     fprintf(2,'\nCopy the following text into the Visual Studio x86 Native Tools Prompt (Start->Programs->MS Visual Studio->Tools):\n\n');
    batstr = '!vcvarsall.bat x86';
    libpath = [path filesep pkgnameA filesep 'Release\'];
else
%     fprintf(2,'\nCopy the following text into the Visual Studio x64 Native Tools Prompt (Start->Programs->MS Visual Studio->Tools):\n\n');
    batstr = '!vcvarsall.bat x86_amd64';
    libpath = [path filesep pkgnameA filesep 'x64\Release\'];
end
outpath = [path filesep pkgnameA filesep 'dump.txt'];

alllibs = [libsA libsB];
str = ['dumpbin /ALL /OUT:"' outpath '" '];
for i = 1:length(alllibs)
    str = [str '"',libpath,alllibs{i} '" '];
end
str = sprintf('%s\n',str);
% %%
cd(vcdir)
eval([batstr ' & ' str]);
cd(cdir);
fprintf('Generated %s containing all library symbols\n',outpath);
%Copy libmwma57 to output dir
if(strcmp(computer,'PCWIN'))
    copyfile([cd '/Solvers/Source/ipopt/lib/libmwma5732.lib'],[libpath filesep 'libmwma57.lib'],'f');
elseif(strcmp(computer,'PCWIN64'))
    copyfile([cd '/Solvers/Source/ipopt/lib/libmwma5764.lib'],[libpath filesep 'libmwma57.lib'],'f');
else
    error('OS not supported!');
end

% Generate .def files to build DLLs
export_def_gen(outpath,[path filesep pkgnameA],pkgnameA);
export_def_gen(outpath,[path filesep pkgnameB],pkgnameB);


%% NOW BUILD optiSolverPackA AND optiSolverPackB In Visual Studio

%Once built, execute this cell to copy outputs
if(strcmp(computer,'PCWIN'))
    copyfile([libpath filesep 'optiSolverPackA32.lib'],[cd '/Solvers/Source/optiSolverPackA32.lib'],'f');
    copyfile([libpath filesep 'optiSolverPackA32.dll'],[cd '/Solvers/optiSolverPackA32.dll'],'f');
    copyfile([libpath filesep 'optiSolverPackB32.lib'],[cd '/Solvers/Source/optiSolverPackB32.lib'],'f');
    copyfile([libpath filesep 'optiSolverPackB32.dll'],[cd '/Solvers/optiSolverPackB32.dll'],'f');
elseif(strcmp(computer,'PCWIN64'))
    copyfile([libpath filesep 'optiSolverPackA64.lib'],[cd '/Solvers/Source/optiSolverPackA64.lib'],'f');
    copyfile([libpath filesep 'optiSolverPackA64.dll'],[cd '/Solvers/optiSolverPackA64.dll'],'f');
    copyfile([libpath filesep 'optiSolverPackB64.lib'],[cd '/Solvers/Source/optiSolverPackB64.lib'],'f');
    copyfile([libpath filesep 'optiSolverPackB64.dll'],[cd '/Solvers/optiSolverPackB64.dll'],'f');
end

%% Build MEX Files

%Files to compile
files = {'clpmex.cpp','coinR.cpp','coinW.cpp','cbcmex.cpp','glpkcc.cpp','mumpsmex.c','mumpsmex.c','ipopt','bonmin',...
         'filtersdmex.c','filtersdmex.c','lbfgsb','levmarmex.c','m1qn3mex.c','hybrjmex.c','lmdermex.c',...
         'nl2solmex.c','ooqpmex.cpp','nloptmex.c','mkltrnls.c','mklJac.c','dsdpmex.c','pswarmmex.c','lpsolve',...
         'nomadmex.cpp'};
%Corresponding output names
outnames = {'clp','coinR','coinW','cbc','glpk','mumps','zmumpsmex','ipopt','bonmin','filtersd','filtersdsp','lbfgsb',...
            'levmar','m1qn3','hybrj','lmder','nl2sol','ooqp','nlopt','mkltrnls','mklJac','dsdp','pswarm','lp_solve',...
            'nomad'};
%Associated Preprocessor Symbols
pp = {'COIN_MSVS','COIN_MSVS','COIN_MSVS','COIN_MSVS',[],'MUMPS_ARITH=2','MUMPS_ARITH=8',...
      {'IPOPT_BUILD','haveMKL','haveMKLPARDISO','haveMA57'},...
      {'BONMIN_BUILD','IPOPT_BUILD','haveMKL','haveMKLPARDISO','haveMA57'},...
      {'STATICLINK','libfilterSD'},{'STATICLINK','libfilterSDsp','SPARSEVER'},[],[],[],[],[],...
      [],{'HAVE_MA57','HAVE_PARDISO'},[],{'MKLLINKSEQ'},{'MKLLINKSEQ'},[],[],{'MATLAB','WIN32','LPSOLVEAPIFROMLIB'},...
      'OPTI_VERSION'};
	  
files = {'nomadmex.cpp'};
outnames = {'nomad'};
pp = {'OPTI_VERSION'};	 

fprintf('\n------------------------------------------------\n');
fprintf('OPTI SOLVER PACK MEX FILE INSTALL\n\n');

%CD to Source Directory
cdir = cd;
cd 'Solvers/Source';

%Setup Common Include Directories 
inc = [' -I' clppath filesep 'src -I' clppath filesep '..\CoinUtils\src -I' glpkpath filesep 'src '];
inc = [inc '-I' cbcpath filesep 'src -I' cbcpath '\..\Cgl\src -I' cbcpath '\..\Osi\src\Osi '];
inc = [inc '-I' clppath filesep 'src\OsiClp -I' mumpspath '\include -I' mumpspath '\PORD\include '];
inc = [inc '-I' ippath '/src/Algorithm -I' ippath '/src/Common -I' ippath '/src/Interfaces '];
inc = [inc '-I' ippath '/src/LinAlg -I' bminpath '/BuildTools/headers '];
inc = [inc '-I' bminpath '/Bonmin/src/Interfaces -I' bminpath '/Bonmin/src/CbcBonmin -I' bminpath '/Bonmin/src/Algorithms '];
inc = [inc '-I' bminpath '/Bonmin/src/Algorithms/OaGenerators -I' bminpath '/Bonmin/src/Interfaces/Ipopt '];
inc = [inc '-I' levmpath ' -I' ooqppath '/src/QpGen -I' ooqppath '/src/Abstract -I' ooqppath '/src/Vector '];
inc = [inc '-I' ooqppath '/src/QpSolvers -I' pswmpath ' -I' nloptpath '/api -I' dsdppath '/Include '];
inc = [inc '-I' lpspath ' -I' nmdpath ' '];
%Get MKL
[~,mkl_for_link,mkl_inc,mkl_lib,mkl_clib] = opti_FindMKL();
inc = [inc '-I"' mkl_inc '" '];
%Get MKL Fortran Libs
lib = [' -L"' mkl_lib '" -L"' mkl_clib '" ' mkl_for_link];
%Get Ctrl-C Library
lib = [lib '-llibut '];

for i = 1:length(files)
    %Remove Existing Solver from Memory
    clear(outnames{i});

    %Set Preprocessors
    post = [' -output ' outnames{i}];
    if(~isempty(pp{i}))
        if(iscell(pp{i}))
            if(strcmp(pp{i}{1},'STATICLINK')) %statically link libs
                post = [post ' -L"' libpath '" -l' pp{i}{2}];
                for j = 3:length(pp{i})
                    post = [post ' -D' pp{i}{j}];
                end
            elseif(strcmp(pp{i}{1},'MKLLINKSEQ')) %non-threaded library
                post = [post opti_FindMKL('seq')];                
            else
                for j = 1:length(pp{i})
                    post = [post ' -D' pp{i}{j}];
                end
            end
        else
            post = [post ' -D' pp{i}];
        end
    end
    
    %Get Import Library
    pk = 'B';
    idx = strfind(libsA,outnames{i}(1:3));
    for j = 1:length(idx)
        if(~isempty(idx{j}))    
            pk = 'A';
            break;
        end
    end
    if(strcmp(computer,'PCWIN'))
        libl = [lib ' -L"' cd '" -loptiSolverPack' pk '32  ']; 
    else
        libl = [lib ' -L"' cd '" -loptiSolverPack' pk '64  '];
    end
        
    %Compile & Move
    if(isempty(strfind(files{i},'.'))) %multiple files
        %Cd to this folder
        cd(files{i});
        %Find all files in here
        f = dir;
        %Create a list of ones to compile
        str = [];
        for j = 1:length(f)
            if(~f(j).isdir && (strcmp(f(j).name(end-2:end),'cpp') || strcmp(f(j).name(end),'c')))
                str = [str ' "' cd filesep f(j).name '" '];
            end
        end
        %Add include folder
        str = [str ' -I"' cd filesep 'Include" '];     
        fprintf('Compiling %s [%s.cpp]...',upper(outnames{i}),files{i});
        pre = ['mex -largeArrayDims ' str];
        cd ../
    else
        fprintf('Compiling %s [%s]...',upper(outnames{i}),files{i});
        pre = ['mex -largeArrayDims ' files{i}];
    end
    try
        eval([pre inc libl post])
        movefile([outnames{i} '.' mexext],'../','f')
        fprintf('Done!\n');
    catch ME
        cd(cdir);
        error('opti:packcompile','Error Compiling %s!\n%s',upper(outnames{i}),ME.message);
    end    
end
%Return to user path
cd(cdir);
fprintf('------------------------------------------------\n');







