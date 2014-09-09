% HYBRJ  Solve a SNLE using HYBRJ
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_hybrj() INSTEAD!
%
% hybrj uses the Minpack Powell Hybrid library.
%
%   [x,fval,exitflag,fevals] = hybrj(fun,grad,x0,opts)
%
%   Input arguments:
%       fun - vector of nonlinear functions (equations to be solved)
%       grad - gradient of nonlinear equations (optional)
%       x0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       fevals - number of function evaluations taken
%
%   Option Fields (all optional):
%       display - solver display level [0,1,2]
%       maxfeval - maximum function evaluations
%       maxtime - maximum solver execution time
%       iterfun - Iteration Callback Function, stop = iterfun(iter,fval,x)
%
%   Return Status:
%       1 - optimal
%       0 - maximum iterations exceeded
%      -1 - infeasible / could not converge
%      -2 - hybrj error / other
%      -5 - user exit
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)