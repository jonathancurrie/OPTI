function [x,fval,exitflag,info] = opti_pswarm(fun,lb,ub,x0,A,b,opts)
%OPTI_PSWARM Solve a Global NLP using PSwarm (Pattern & Particle Swarm)
%
%   min f(x)       subject to:   A*x <= b
%    x                           lb <= x <= ub
%
%   x = opti_pswarm(fun,lb,ub,x0) solves a Global NLP where fun is the 
%   objective function and lb and ub are the decision variable bounds 
%   (required and must be finite). x0 is a starting guess.
%
%   x = opti_pswarm(fun,lb,ub,x0,A,b) solves subject to the linear
%   inequality constraints A*x <= b.
%
%   x = opti_pswarm(fun,...,b,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_pswarm(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR PSwarm
%   See referenced Lesser GNU Public License

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(nargin < 7), opts = optiset; else opts = optiset(opts); end
if(nargin < 6), b = []; end
if(nargin < 5), A = []; end
if(nargin < 4), error('PSwarm requires at least 4 arguments'); end

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('PSWARM requires an initial guess, x0!');
end

%Add in pswarmset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = pswarmset(opts.solverOpts);
    %Add in options from optiset
    popts.maxiter = opts.maxiter;
    popts.maxfeval = opts.maxfeval;
    popts.maxtime = opts.maxtime;
    popts.tolrfun = opts.tolrfun;
    popts.iterfun = opts.iterfun;
else
    popts = opts;
end
%Setup display level
popts.display = dispLevel(opts.display);

%Ensure we have bounds
if(isempty(lb) || isempty(ub))
    error('PSwarm only solves bounded NLPs (or bounded NLS/SCNLE) - You must supply lb and ub to this function (or the OPTI constructor).');
end
%Ensure we have finite bounds
if(any(isinf(lb)) || any(isinf(ub)))
    error('PSwarm requires all bounds to be finite');
end

t = tic;
% Run PSwarm
[x, fval, exitflag, iter, feval] = pswarm(fun,x0,lb,ub,A,b,popts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'PSwarm: Pattern & Particle Swarm Global Optimization';

switch(exitflag)
    case 1
        info.Status = 'Converged';
    case 0
        info.Status = 'Exceeded Iterations / Function Evaluations / Time';
    case -1
        info.Status = 'Abnormal Exit';
    case -2
        info.Status = 'PSwarm Error';
    case -5
        info.Status = 'User Exit';
    otherwise        
        info.Status = 'PSwarm Error';
end
