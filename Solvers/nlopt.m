% NLOPT  Solve a NLP using NLOPT
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_nlopt() INSTEAD!
%
% nlopt uses the NLOPT library (LGPL).
%
%   [x,fval,status] = nlopt(nlprob,x0)
%
%   Input arguments:
%       nlprob - problem description structure (see below)
%       x0 - initial solution guess
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       status - exit status (see below)
%       funcevals - structure of function and constraint #evals
%
%   nlprob Fields:
%       objective - objective function handle
%       gradient - objective gradient function handle (optional, req for *D algorithms)
%       nlcon - nonlinear constraints function handle (both ineq and eq)
%       nlrhs - nonlinear constraints rhs (nlcon(x) <=,==,>= nlrhs)
%       nle - nonlinear constraint type (-1 <=, 0 ==, 1 >=)
%       nljac - nonlinear constraint jacobian function handle (optional, req for *D algorithms)
%       lb - decision variables lower bounds
%       ub - decision variables upper bounds
%       options - options structure (see below)
%
%   options Fields:
%       algorithm - NLOPT algorithm ID [required, integer]
%       maxtime - maximum solver execution time
%       maxfeval - maximum solver evaluation time
%       tolrfun - relative function tolerance
%       tolafun - absolute function tolerance
%       display - 0 to 2 display level
%       iterfun - callback function for each function evaluation [iterfun(nfeval,f,x)]
%       local_optimizer - structure containing settings for local
%       optimizer, includes fields algorithm, maxtime, maxfeval, tolrfun
%       and tolafun as above.
%
%   Return Status:
%       1 - generic success
%       2 - success (stopval reached)
%       3 - success (ftol_rel or ftol_abs reached)
%       4 - success (xtol_rel or xtol_abs reached)
%       5 - max fevals reached
%       6 - max time reached
%       -1 - generic failure
%       -2 - invalid arguments (check for bounds on global solvers)
%       -3 - out of memory
%       -4 - halted due to roundoff errors
%       -5 - user exit
%
 
%   Based in parts on the original MEX interface by Prof. Steven Johnson.
%
%   Copyright (C) 2013 Jonathan Currie (IPL)