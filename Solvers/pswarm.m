% PSWARM  Solve a Global NLP using PSWARM
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_pswarm() INSTEAD!
%
% pswarm uses the Pattern and Particle Swarm Optimization library.
%
%   [x,fval,exitflag,iter,feval] = pswarm(fun,x0,lb,ub,A,b,opts)
%
%   Input arguments:
%       fun - nonlinear function handle
%       x0 - initial solution guess
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       A - linear inequality matrix (dense)
%       b - linear inequality rhs
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       feval - number of function evaluations taken by the solver
%
%   Option Fields (all optional - see pswarmset):
%       display - solver display level [0,1,2]
%       tolfun - function tolerance
%       maxiter - maximum solver iterations
%       maxfeval - Maximum number of function evaluations 
%       maxtime - Maximum solver execution time
%       swarm_size - Swarm Size
%       vectorized - Objective function is vectorized
%       mu - Cognitial Parameter 
%       nu - Social Parameter
%       iweight - Initial Weight
%       fweight - Final Weight
%       delta - Initial Delta
%       idelta - Increase Delta
%       ddelta - Decrease Delta
%
%   Return Status:
%       1 - converged
%       0 - maximum iterations / function evaluations exceeded
%      -1 - abnormal exit
%      -2 - memory error
%      -3 - population error
%
%   
%   Copyright (C) 2012 Jonathan Currie (IPL)