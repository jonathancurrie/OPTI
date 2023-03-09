% M1QN3  Solve a UNO using M1QN3
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_m1qn3() INSTEAD!
%
%   [x,fval,exitflag,iter,feval] = m1qn3(fun,grad,x0,opts)
%
%   Input arguments:
%       fun - nonlinear function handle
%       grad - nonlinear gradient function handle (required)
%       x0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       feval - number of function evaluations taken by the solver
%
%   Option Fields (all optional):
%       display - solver display level [0,1,2]
%       tolafun - gradient convergence tolerance
%       maxiter - maximum solver iterations {1000}
%       maxfeval - Maximum number of function evaluations {1500}
%       maxtime - Maximum solver execution time {1000s}
%       nupdate - Number of L-BFGS updates for Hessian (5-10) {5}
%       iterfun - Iteration callback function handle
%
%   Return Status:
%       0 - user exited
%       1 - successful termination
%       2 - input argument error
%       3 - line-search is blocked
%       4 - maximum iterations reached
%       5 - maximum function evaluations reached
%       6 - stop on dxmin during the line-search
%       7 - numerical error
%       8 - maximum time reached
%
%   
%   Copyright (C) 2012 Jonathan Currie (Control Engineering)