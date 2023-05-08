function [x,fval,exitflag,info] = opti_gsl_nls(fun,grad,x0,ydata,opts)
%OPTI_GSL_NLS Solve a NLS using GSL Multifit Nonlinear 
%
%   min sum[ (F(x) - ydata)^2 ] 
%    x
%
%   x = opti_gsl_nls(fun,grad,x0,ydata) solves a NLS where fun is the 
%   fitting function. grad is an optional gradient of the fitting function 
%   and x0 is a starting guess. ydata is the data to fit the function to. 
%
%   x = opti_gsl_nls(fun,...,ydata,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_gsl_nls(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR gsl_multifit_nlinear

%   Copyright (C) 2017 Jonathan Currie (Control Engineering)

if(nargin < 5), opts = optiset; end
if(nargin < 4), error('OPTI_GSL_NLS requires at least 4 arguments'); end

% Setup display level
opts.display = dispLevel(opts.display);
opts.optiver = optiver;

% Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('OPTI_GSL_NLS requires an initial guess, x0!');
end

% Addin gslset settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    sopts = gslset(opts.solverOpts);    
else
    sopts = [];
end
% Add OPTI Options
sopts.maxiter   = opts.maxiter;
sopts.maxfeval  = opts.maxfeval;
sopts.maxtime   = opts.maxtime;
sopts.display   = opts.display;
sopts.tolafun   = opts.tolafun;
sopts.iterfun   = opts.iterfun;

% Construct problem structure
nlprob.fun      = fun;
nlprob.grad     = grad;
nlprob.ydata    = ydata;
nlprob.x0       = x0;
nlprob.options  = sopts;
nlprob.probType = 'nls';

t = tic;
% Run GSL
[x, fval, exitflag, stats] = gsl(nlprob);

%Collect Results
info.Iterations = stats.niter;
info.FuncEvals = stats.nfeval;
info.GradEvals = stats.ngeval;
info.Time = toc(t);
info.Covar = stats.covar;
info.Algorithm = stats.algorithm;

switch(exitflag)
    case 0
        info.Status = 'Success';
        exitflag    = 1;
    case 27
        info.Status = 'No Further Progress Could Be Made';
        exitflag    = -1;
    case 11
        info.Status = 'Exceeded Maximum Iterations';
        exitflag    = 0;
    case -6
        info.Status = 'Exceeded Maximum Time';
        exitflag    = 0;
    case -7
        info.Status = 'Exceeded Maximum Function Evaluations';
        exitflag    = 0;
    case -8
        info.Status = 'No Progress in First Iteration';
        exitflag    = -3;
    case -5
        info.Status = 'User Exited';
        exitflag    = -5;
    otherwise        
        info.Status = sprintf('GSL Error (Code %d)', exitflag);
        exitflag    = -2;
end
