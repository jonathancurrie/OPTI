function [x,fval,exitflag,info] = opti_hybrj(fun,grad,x0,opts)
%OPTI_HYBRJ Solve a SNLE using HYBRJ (Powell-Hybrid MINPACK Routine)
%
%   F(x) = 0
%
%   x = opti_hybrj(fun,grad,x0) solves a SNLE where fun is a vector of
%   functions to solve. grad is an optional gradient of the nonlinear
%   equations and x0 is a starting guess.
%
%   x = opti_hybrj(fun,grad,x0,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_hybrj(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR HYBRJ + HYBRD
%   See supplied License

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(nargin < 4), opts = optiset; end
if(nargin < 3), error('HYBRJ requires at least 3 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('HYBRJ requires an initial guess, x0!');
end

t = tic;
% Run HYBRJ
[x, fval, exitflag, iter] = hybrj(fun,grad,x0,opts);

%Collect Results
info.FuncEvals = iter;
info.Time = toc(t);
info.Algorithm = 'HYBRJ: MINPACK Powell Hybrid';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Function Evaluations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'HYBRJ Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'HYBRJ Error';
end
