% IPOPT  Solve a NLP using IPOPT
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_ipopt() INSTEAD!
%
% ipopt uses the Interior-Point Optimization library.
%
%   [x,output] = ipopt(x0,funcs,opts)
%
%   Input arguments:
%       x0 - initial solution guess
%       funcs - functions structure (see below)
%       opts - options structure (see below)
%
%   Return arguments:
%       x - solution vector
%       output - return information structure
%
%   Funcs Fields:
%       objective - objective function handle
%       gradient - objective gradient function handle
%       constraints - constraints function handle (optional)
%       jacobian - constraints Jacobian function handle (sparse)
%       jacobianstructure - constraints Jacobian sparsity structure (sparse)
%       hessian - objective Hessian of the Lagrangian (sparse,optional)
%       hessianstructure - objective Hessian sparsity structure (sparse)
%       
%   Opts Fields:
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       cl - constraints lower bounds
%       cu - constraints upper bounds
%       ipopt - solver settings structure (see ipoptset.m)
%
%   Return Status:
%       {0,1} - optimal
%      -1 - exceeded iterations
%       2 - infeasible
%       {3,4,5,-2,-3,-10} - unbounded or infeasible
%       other - ipopt error / other
%
%   
%   Based almost entirely on the original MEX interface by Dr. Peter Carbonetto.
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)