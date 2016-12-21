function [x,fval,exitflag,info,Opt] = opti_fminunc(fun,x0,opts)
%OPTI_FMINUNC Solve an UNO using an OPTI UNO Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_fminunc(fun,x0) solves the unconstrained 
%   nonlinear optimization min f(x) where fun is the nonlinear function to 
%   be minimized [fun(x)], starting at x0.
%
%   [x,fval,exitflag,info] = opti_fminunc(fun,...,opts) allows the 
%   user to specify optiset options. This includes specifying a solver via 
%   the 'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_fminunc(fun,...) returns the internally 
%   built OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 3, opts = optiset; end 
if nargin < 2, error('You must supply at least 2 arguments to opti_fminunc'); end

%Sort out Fun + Grad
[f,g] = detGrad(fun,x0);

%Build OPTI Object
Opt = opti('fun',f,'grad',g,'x0',x0,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);
