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
