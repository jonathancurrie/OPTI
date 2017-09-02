function [x,fval,exitflag,info] = opti_gsl(nlprob, x0)
%OPTI_GSL Solve a NLS using GSL (GNU Scientific Library)
%
%   x = opti_gsl(nlprob, x0) solves a...
%
%   [x,fval,exitflag,info] = opti_gsl(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR GSL
%   See supplied License

%   Copyright (C) 2017 Jonathan Currie (IPL)

if (nargin < 2), error('GSL requires at least 2 arguments'); end
if (~isfield(nlprob,'options') || isempty(nlprob.options))
    nlprob.options = optiset;
end

%Setup display level
nlprob.options.display = dispLevel(nlprob.options.display);
nlprob.options.optiver = optiver;

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('GSL requires an initial guess, x0!');
end
nlprob.x0 = x0;

t = tic;
% Run GSL
[x, fval, exitflag, stats] = gsl(nlprob);

%Collect Results
info.Iterations = stats.niter;
info.FuncEvals = stats.nfeval;
info.GradEvals = stats.ngeval;
info.Time = toc(t);
% info.Algorithm = ['GSL: ' stats.algorithm];
info.Covar = stats.covar;

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Function Evaluations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'GSL Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'GSL Error';
end
