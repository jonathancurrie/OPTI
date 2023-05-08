% GSL  Solve a NLS using GSL
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_gsl_xxx() INSTEAD!
%
% gsl uses the GNU Scientific Library.
%
%   [x,fval,exitflag,iter] = gsl(prob)
%
%   Input arguments:
%       prob - a structure containing the problem definition (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%
%   Problem Structure Fields (problem type specific):
%       NLS:
%           probType    - must be specified as 'nls'
%           fun         - nonlinear fitting function handle
%           grad        - gradient of nonlinear fitting function handle (optional)
%           x0          - initial solution guess
%           ydata       - fitting data
%           opts        - solver options (see below)
%
%   Option Fields (all optional - see gslset for problem specific):
%       display - solver display level [0,1,2]
%       maxiter - maximum solver iterations
%       maxfeval - maximum function evaluations
%       maxtime - maximum solver execution time
%       tolafun - absolute function tolerance
%       iterfun - Iteration Callback Function, stop = iterfun(iter,fval,x)
%
%   Return Status:
%       0 - success 
%      11 - maximum iterations exceeded
%      -6 - maximum time exceeded
%      -7 - maximum function evaluations exceeded
%      27 - no progress could be made
%      -8 - no progress in the first iteration
%      -5 - user exited
%   
%
%   Copyright (C) 2017 Jonathan Currie (Control Engineering)