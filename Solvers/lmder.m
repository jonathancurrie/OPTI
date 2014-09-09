% LMDER  Solve a NLS using LMDER
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_lmder() INSTEAD!
%
% lmder uses the Minpack Levenberg-Marquardt library.
%
%   [x,fval,exitflag,feval] = lmder(fun,grad,x0,ydata,opts)
%
%   Input arguments:
%       fun - nonlinear fitting function handle
%       grad - gradient of nonlinear fitting function handle (optional)
%       x0 - initial solution guess
%       ydata - fitting data
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       feval - number of function evaluations
%
%   Option Fields (all optional):
%       display - solver display level [0,1,2]
%       tolrfun - function tolerance
%       maxfeval - maximum function evaluations
%       maxtime - maximum solver execution time       
%       iterfun - Iteration Callback Function, stop = iterfun(iter,fval,x)
%
%   Return Status:
%       1 - optimal
%       0 - maximum iterations exceeded
%      -1 - infeasible / could not converge
%      -2 - lmder error / other
%      -5 - user exit
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)