%% OPTI Toolbox File IO Examples
%
% This file loads a number of supplied examples and shows how to solve
% them using the OPTI Toolbox. You should read and complete Basic_Usage.m 
% BEFORE running the below examples.
%
% The underlying parser uses CoinUtils + GLPK for reading and writing all
% file types.
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

%This page is supplemented by examples on the following pages:
web('https://www.inverseproblem.co.nz/OPTI/index.php/File/MPS');
web('https://www.inverseproblem.co.nz/OPTI/index.php/File/SDPA','-new');
web('https://www.inverseproblem.co.nz/OPTI/index.php/File/GAMS','-new');

%% Loading a MPS Problem
% OPTI Toolbox is supplied with a number of example LP, MILP & QP problems
% in MPS, QPS, MOD and LP format, located in Test Problems/. To
% load a MPS problem, simply use the command below. Returned will be a 
% optiprob structure containing the data in the MPS file. Note all matrices 
% are returned as sparse matrices. If you do not specify a file extension 
% you must supply the  file type via a second argument.
clc
prob = coinRead('testLP.mps')

%% Example 1 - Solving a Loaded MPS Problem
% Solving a loaded MPS problem is simple, just pass it to the opti
% constructor and call solve:
clc
Opt = opti(prob) 

x = solve(Opt)

%% Loading a SDPA Model
% You can also load a problem by specifying the file type explicitly, as
% well as printing the headers + errors as they exist:
clc
prob = sdpRead('arch0','sdpa-s')

%% Example 2 - Solving a Loaded QPS Problem
% Solving a loaded QPS problem is the same as an MPS problem:
clc
Opt = opti(prob,optiset('solver','csdp','maxiter',5e3)) 

[x,fval,e,info] = solve(Opt);
info


%% Example 5 - Writing a MPS/QPS Problem
% You can also write an opti problem to most file types. The following
% example writes the an AMPL model to a GAMS model

prob = amplRead('diet.nl');
Opt = opti(prob)

write(Opt,'testGAMSLP.gms')

%% NOTE - rememeber to delete the output files!
delete testgamslp.gms

