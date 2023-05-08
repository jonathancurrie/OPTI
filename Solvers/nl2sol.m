% NL2SOL  Solve a NLS using NL2SOL
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_nl2sol() INSTEAD!
%
% nl2sol uses the Adaptive Nonlinear Least Squares library.
%
%   [x,fval,exitflag,iter,feval] = nl2sol(fun,grad,x0,ydata,lb,ub,opts)
%
%   Input arguments:
%       fun - nonlinear fitting function handle
%       grad - gradient of nonlinear fitting function handle (optional)
%       x0 - initial solution guess
%       ydata - fitting data
%       lb - decision variable lower bounds (optional)
%       ub - decision variable upper bounds (required if lb present)
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       feval - number of function evaluations
%
%   Option Fields (all optional):
%       display - solver display level [0,1,2]
%       tolrfun - relative function tolerance
%       tolafun - absolute function tolerance
%       maxiter - maximum solver iterations
%       maxfeval - maximum solver function evaluations
%
%   Return Status:
%       1 - optimal / stopped by xtol or ftol or both
%       0 - maximum iterations exceeded
%      -1 - could not converge / tolerance too small
%      -2 - nl2sol error / bad parameters / out of range
%      -3 - unknown error
%      -5 - user exit
%
%
%   Copyright (C) 2011 Jonathan Currie (Control Engineering)