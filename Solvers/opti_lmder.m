function [x,fval,exitflag,info] = opti_lmder(fun,grad,x0,ydata,opts)
%OPTI_LMDER Solve a NLS using LM_DER (Levenberg-Marquardt MINPACK Routine)
%
%   min sum[ (F(x) - ydata)^2 ] 
%    x
%
%   x = opti_lmder(fun,grad,x0,ydata) solves a NLS where fun is the fitting
%   function. grad is an optional gradient of the fitting function and x0 
%   is a starting guess. ydata is the data to fit the function to. Use
%   ydata = zeros() for solving a SNLE.
%
%   x = opti_lmder(fun,grad,x0,ydata,opts) uses opts to pass optiset 
%   options to the solver. 
%
%   [x,fval,exitflag,info] = opti_lmder(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR LM_DER + LM_DIF
%   See supplied License

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(nargin < 5), opts = optiset; end
if(nargin < 4), error('LM_DER requires at least 4 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('LMDER requires an initial guess, x0!');
end

t = tic;
% Run LM_DER
[x, fval, exitflag, feval] = lmder(fun,grad,x0,ydata,opts);

%Collect Results
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'LM_DER: MINPACK Levenberg-Marquardt';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Function Evaluations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'LM_DER Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'LM_DER Error';
end
