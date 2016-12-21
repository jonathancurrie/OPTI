% FILTERSDSP  Solve a SPARSE NLP using FilterSD
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_filtersd() INSTEAD!
%
% FilterSD is a NLP solver written in FORTRAN and available from:
% https://projects.coin-or.org/filterSD
%
%   [x,fval,exitflag,stats,lambda] = filtersdsp(fun, grad, x0, lb, ub, nlcon, nljac, nljacstr, cl, cu, opts)
%
%   Input arguments:
%       fun - nonlinear function handle
%       grad - gradient of nonlinear function handle
%       x0 - initial solution guess
%       lb - decision variable lower bounds 
%       ub - decision variable upper bounds 
%       nlcon - nonlinear constraints handle 
%       nljac - Jacobian of nonlinear constraints handle (sparse matrix)
%       nljacstr - Sparse matrix containing locations of ALL possible non-zero locations in Jacobian
%       cl - nonlinear constraints lower bound
%       cu - nonlinear constraints upper bound
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       stats - statistics structure
%       lambda - Lagrange multipliers at solution [bounds;nl constraints]
%
%   Option Fields (all optional):
%       maxiter - maximum major solver iterations
%       maxtime - maximum solver execution time [won't return last solution on bounded problems]
%       maxfeval - maximum function evaluations [won't return last solution on bounded problems]
%       display - solver display level [0-2]
%       rho - initial trust region radius
%       htol - tolerance allowed in sum h(x) of constraint infeasibilities
%       rgtol - tolerance allowed in reduced gradient L2 norm
%       iterfun - iteration callback of the form ifun(feval,fval,x)
%
%   Return Status:
%       0 - Successful run
%       1 - Unbounded NLP (f(x) <= fmin at htol feasible point)
%       2 - Bounds on x are inconsistent
%       3 - Local minimum of feasibility problem and h(x) > htol (nonlinear constraints are locally inconsistent)
%       4 - Initial point x has h(x) > ubd
%       5 - Iteration limit reached
%       6 - Termination with rho <= htol
%       7 - Not enough workspace memory in ws or lws
%       8 - Insufficient space for filter
%       9 and above - LCP solver problem
%       105 - User Exit
%
%   Copyright (C) 2013 Jonathan Currie (IPL)