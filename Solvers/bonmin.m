% BONMIN  Solve a MINLP using BONMIN
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_bonmin() INSTEAD!
%
% bonmin uses the Basic Open source Nonlinear Mixed Integer library.
%
%   [x,output] = bonmin(x0,funcs,opts)
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
%       ipopt - relaxed solver settings structure (see ipoptset.m)
%       bonmin - bonmin settings (see bonminset.m)
%       var_type - integer variable indicies (-1 binary, 0 cont, 1 integer)
%       var_lin - decision variable linearity (0 nonlinear, 1 linear)
%       cons_lin - constraint linearity (0 nonlinear, 1 linear)
%       display - log level (0 off, 1-5 increasing)
%
%   Return Status:
%       0 - proven optimal
%      -1 - exceeded iterations
%       1 - proven infeasible
%       2 - feasible
%       3 - unbounded or infeasible
%       4 - no solution known
%       other - bonmin error / other
%
%   
%   Based almost entirely on the MEX interface to IPOPT by Dr. Peter Carbonetto.
%
%   Copyright (C) 2011 Jonathan Currie (IPL)