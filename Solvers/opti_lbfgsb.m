function [x,fval,exitflag,info] = opti_lbfgsb(fun,grad,lb,ub,x0,opts)
%OPTI_LBFGSB Solve a bounded NLP using L-BFGS-B (Nocedal & Zhu)
%
%   min f(x)       subject to:   lb <= x <= ub
%    x
%
%   x = opti_lbfgsb(fun,grad,lb,ub,x0) solves a NLP where fun is the 
%   objective function and grad is the gradient of the objective. lb and ub
%   are the decision variable bounds (required) and x0 is a starting guess.
%   Avoid Infinite bounds.
%
%   x = opti_lbfgsb(fun,...,x0,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_lbfgsb(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR L-BFGS-B
%   See referenced BSD License

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(nargin < 6), opts = optiset; end
if(nargin < 5), error('LBFGSB requires at least 5 arguments'); end

%Add in clpset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = lbfgsbset(opts.solverOpts);    
else
    popts = lbfgsbset;
end
%Add in options from optiset   
popts.maxiter = opts.maxiter;
popts.maxtime = opts.maxtime;
popts.maxfeval = opts.maxfeval;
popts.tolrfun = opts.tolrfun;
popts.iterfun = opts.iterfun;
popts.optiver = optiver;
%Setup display level
popts.display = dispLevel(opts.display);

%Ensure we have bounds
if(isempty(lb) || isempty(ub))
    error('L-BFGS-B only solves bounded NLPs - You must supply lb and ub to this function');
end
%Ensure we have a gradient
if(isempty(grad))
    error('L-BFGS-B requires a gradient function');
end

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('L-BFGS-B requires an initial guess, x0!');
end

t = tic;
% Run L-BFGS-B
[x, fval, status, iter, feval] = lbfgsb(fun,grad,lb,ub,x0,popts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'L-BFGS-B: Limited Memory BFGS Bounded Optimization';

switch(status)
    case 0
        info.Status = 'Success';
        exitflag = 1;
    case 1
        info.Status = 'Abnormal Termination';
        exitflag = -2;
    case 2
        info.Status = 'Error On Input';
        exitflag = -3;
    case 3
        info.Status = 'Exceeded Iterations';
        exitflag = 0;
    case 4
        info.Status = 'Exceeded Maximum Fevals';
        exitflag = 0;
    case 5
        info.Status = 'Exceeded Maximum Time';
        exitflag = 0;
    case 6
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise        
        info.Status = 'L-BFGS-B Error';
        exitflag = -3;
end
