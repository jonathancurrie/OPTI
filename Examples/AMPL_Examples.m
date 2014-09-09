%% OPTI Toolbox AMPL Examples
%
% This file loads a number of supplied examples and shows how to solve
% them using the OPTI Toolbox. You should read and complete BasicUsage.m 
% BEFORE running the below examples.
%
% The underlying parser uses Netlib's AMPL Solver Library (ASL) Interface.
%
%   Copyright (C) 2014 Jonathan Currie (I2C2)

% There is also a page on the Wiki which supplements this example:
web('http://i2c2.aut.ac.nz/Wiki/OPTI/index.php/File/AMPL');

%% Loading an AMPL Problem
% OPTI Toolbox is supplied with a number of example AMPL problems ranging
% from LP to NLP including continuous and discrete problems. To load an
% AMPL model it must be in .NL format, see the supplied user's guide page
% on AMPL Interfacing in order to generate the .NL file from your model. To
% load an AMPL problem, simply use the command below. Returned will be a 
% optiprob structure containing the data in the file. 
clc
prob = amplRead('diet.nl')

%% Example 1 - Solving a Loaded AMPL Problem
% Solving a loaded AMPL problem is simple, just pass it to the opti
% constructor and call solve. You will note that the command amplRead will
% automatically determine the problem type (LP, QP, NLP, etc) based on
% properties supplied by the ASL.
clc
Opt = opti(prob) 

x = solve(Opt)

%% Loading an AMPL Nonlinear Program
% The AMPL format can easily specify nonlinear problems unlike the LP or
% MPS formats. 
clc
prob = amplRead('ch3.nl')

%% Example 2 - Solving a Nonlinear Problem from AMPL
% In order to solve a nonlinear problem the model .nl file
% remains open after the amplRead function is called to allow callbacks for
% the objective, gradient, etc. Calling solve() will automatically close
% the interface for linear and quadratic problems, however you must manually
% call close for nonlinear problems. Alternatively it will be automatically 
% called once Matlab is closed to clean up ASL memory.
clc
Opt = opti(prob)

x = solve(Opt)

asl('close')

%% Loading an Integer Problem
% AMPL will reorder the variables from the original model based on internal
% rules (see the supplied pdfs). Based on this the ASL indicates which
% variables have binary or integer constraints, which are passed through
% the amplRead interface:
clc
prob = amplRead('multmip1.nl')

%% Example 3 - Solving an AMPL Integer Problem
% As the optiprob structure can fully define an integer problem, no change
% is required to build and solve the MIP
clc
Opt = opti(prob) 

[x,fval] = solve(Opt);
fval

%% A Note On Hessians
% As of OPTI v1.58 the Hessian of an NLP will automatically be added to the
% problem structure. However note it is the Hessian of the Lagrangian, with
% the calling form:
%
%   H = hessian(x,sigma,lambda)
%
% This means the Hessian includes second derivatives of both the objective
% and the constraints. Sigma is a scalar scaling factor on the objective
% derivatives only, and lambda are the lagrange multipliers at the current
% point x (scaling the constraint derivatives).
