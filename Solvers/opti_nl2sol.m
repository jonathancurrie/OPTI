function [x,fval,exitflag,info] = opti_nl2sol(fun,grad,x0,ydata,lb,ub,opts)
%OPTI_NL2SOL Solve a NLS using NL2SOL (Adaptive Nonlinear Least Squares)
%
%   min sum[ (F(x) - ydata)^2 ]       subject to:   lb <= x <= ub
%    x
%
%   x = opti_nl2sol(fun,grad,x0,ydata) solves a NLS where fun is the fitting
%   function. grad is an optional gradient of the fitting function and x0 
%   is a starting guess. ydata is the data to fit the function to. Use
%   ydata = zeros() for solving a SNLE.
%
%   x = opti_nl2sol(fun,grad,x0,ydata,lb,ub) solves subject to decision
%   variables bounds lb <= x <= ub. Infinite bounds are OK.
%
%   x = opti_nl2sol(fun,grad,x0,ydata,lb,ub,opts) uses opts to pass optiset 
%   options to the solver. 
%
%   [x,fval,exitflag,info] = opti_nl2sol(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR NL2SOL + NL2SNO
%   See referenced ACM License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(nargin < 7), opts = optiset; end
if(nargin < 6), ub = []; end
if(nargin < 5), lb = []; end
if(nargin < 4), error('NL2SOL requires at least 4 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('NL2SOL requires an initial guess, x0!');
end

t = tic;
% Run NL2SOL
[x, fval, exitflag, iter, feval] = nl2sol(fun,grad,x0,ydata,lb,ub,opts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'NL2SOL: Adaptive Nonlinear Least Squares';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'NL2SOL Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'NL2SOL Error';
end
