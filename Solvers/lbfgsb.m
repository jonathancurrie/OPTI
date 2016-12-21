% LBFGSB  Solve a Bounded NLP using LBFGSB
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_lbfgsb() INSTEAD!
%
% lbfgsb uses the Limited Memory BFGS Bounded Optimization library.
%
%   [x,fval,exitflag,iter] = lbfgsb(fun,grad,lb,ub,x0,opts)
%
%   Input arguments:
%       fun - nonlinear function handle
%       grad - gradient of nonlinear function handle
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       x0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%
%   Option Fields (all optional):
%       tolrfun - relative function tolerance {1e-7}
%       pgtol   - projected gradient tolerance {1e-5}
%       nupdate - number of L-BFGS updates for Hessian (5-20) {5}
%       maxiter - maximum iterations {1000}
%       maxfeval - maximum function (and gradient) evaluations {1500}
%       maxtime - maximum solver execution time {1000s}
%       display - solver display level [0,1,2]
%       iterfun - iteration callback function, stop = iterfun(iter,fval,x)
%
%   Return Status:
%       0 - success
%       1 - abnormal termination
%       2 - error on input
%       3 - exceeded iterations
%       4 - exceeded max fevals
%       5 - exceeded max time
%       6 - user exit
%   
%   Based in parts on the original MEX interface by Dr. Peter Carbonetto.
%
%   Copyright (C) 2013 Jonathan Currie (IPL)