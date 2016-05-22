%% Miscellaneous Solver/Utility Install for OPTI Toolbox

%   Copyright (C) 2015 Jonathan Currie (I2C2)

% This file will help you compile a few miscellaneous solvers and
% utilities for use with OPTI.

% My build platform:
% - Windows 7 x64
% - Visual Studio 2015

%% NETLIB BLAS
% Download from http://www.netlib.org/blas/

% Build VS Ifort Solution to build .lib
path = 'C:\Solvers\blas-3.6.0';
opti_VSBuild('blas',path)


%% NETLIB LAPACK
% Download from http://www.netlib.org/lapack/

% Build VS Ifort Solution to build .lib
path = 'C:\Solvers\lapack-3.6.0';
opti_VSBuild('lapack',path)


%% MA27 (HSL License - Do not redistribute)
% MA27 is an older sparse linear solver. Download from HSL:
% http://www.hsl.rl.ac.uk/download/MA27/1.0.0/a/

% Build VS Ifort Solution to build .lib [32bit int version]
path = 'C:\Solvers\ma27';
opti_VSBuild('ma27',path)


%% MA57 (HSL License - Do not redistribute)
% MA57 is a new sparse linear solver which replaces MA27. Download from HSL:
% http://www.hsl.rl.ac.uk/download/MA57/3.9.0/

% Build VS Ifort Solution to build .lib [32bit int version]
path = 'C:\Solvers\ma57-3.7.0';
metispath = 'C:\Solvers\metis-4.0.3'; % FULL path to METIS [max version 4.0.3]
opts = []; opts.expaths = metispath;
opti_VSBuild('ma57',path,opts)


%% MATLAB's MA57 Wrapper
% Source is included with OPTI, just call opti_VSBuild as below

opti_VSBuild('mwma57')