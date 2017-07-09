function [x,fval,exitflag,info] = opti_mkltrnls(fun,grad,x0,ydata,lb,ub,opts)
%OPTI_MKLTRNLS Solve a NLS using MKLTRNLS (Intel MKL Trust Region Solver)
%
%   min sum[ (F(x) - ydata)^2 ]       subject to:   lb <= x <= ub
%    x
%
%   x = opti_mkltrnls(fun,grad,x0,ydata) solves a NLS where fun is the 
%   fitting function. grad is an optional gradient of the fitting function 
%   and x0 is a starting guess. ydata is the data to fit the function to. 
%
%   x = opti_mkltrnls(fun,grad,x0,ydata,lb,ub) solves subject to decision
%   variables bounds lb <= x <= ub. Avoid Infinite bounds.
%
%   x = opti_mkltrnls(fun,...,ub,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_mkltrnls(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR DTRNLS, DTRNLSBC and DJACOBI

%   Copyright (C) 2011 Jonathan Currie (IPL)

if(nargin < 7), opts = optiset; end
if(nargin < 6), ub = []; end
if(nargin < 5), lb = []; end
if(nargin < 4), error('MKLTRNLS requires at least 4 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);
opts.optiver = optiver;

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('MKLTRNLS requires an initial guess, x0!');
end

%Addin mkltrnls settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    mopts = mkltrnlsset(opts.solverOpts);    
else    
    mopts = [];
end
%Add OPTI Options
mopts.maxiter = opts.maxiter;
mopts.maxtime = opts.maxtime;
mopts.display = dispLevel(opts.display);
mopts.tolrfun = opts.tolrfun;
mopts.tolafun = opts.tolafun;
mopts.iterfun = opts.iterfun;

t = tic;
%MKLTRNLS requires x0 within bounds
x0 = movex0(lb,ub,x0);
% Run MKLTRNLS
[x, fval, exitflag, iter, feval, geval] = mkltrnls(fun,grad,x0,ydata,lb,ub,mopts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.GradEvals = geval;
info.Time = toc(t);
info.Algorithm = 'MKLTRNLS: Intel MKL Trust Region NLS';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'Singular / Error';
    case -5
        info.Status = 'User Exit';
    otherwise        
        info.Status = 'MKLTRNLS Error';
end
