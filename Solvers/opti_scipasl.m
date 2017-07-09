function [x,fval,exitflag,info] = opti_scipasl(file,opts)
%OPTI_SCIPASL Solve a NLP/MINLP using SCIP to Global Optimality using the 
%             AMPL Interface
%
%   min fun(x)                 subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                              cl <= nlcon(x) <= cu
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1} [i!=j]
%
%   Full Calling Form:
%     [x,fval,exitflag,info] = opti_scipasl(file,opts)
%
%   x = opti_scipasl(file) solves the AMPL model supplied as 'file'. This
%   must be an AMPL .nl model.
%
%   x = opti_scipasl(file,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_scipasl(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR SCIP USING THE MEX INTERFACE
%   See supplied ZIB Academic License

%   Copyright (C) 2012/2013 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 2, opts = optiset('warnings','off'); end
if nargin < 1, error('You must supply at least one argument to opti_scipasl'); end

%Addin scip settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    sopts = scipset(opts.solverOpts);    
else    
    sopts = [];
end
%Add OPTI Options
if(isfield(opts,'maxtime') && ~isempty(opts.maxtime))
    sopts.maxtime = opts.maxtime;
end
if(isfield(opts,'maxiter') && ~isempty(opts.maxiter))
    sopts.maxiter = opts.maxiter;
end
if(isfield(opts,'maxnodes') && ~isempty(opts.maxnodes))
    sopts.maxnodes = opts.maxnodes;
end
if(isfield(opts,'tolrfun') && ~isempty(opts.tolrfun))
    sopts.tolrfun = opts.tolrfun;
end
if(isfield(opts,'objbias') && ~isempty(opts.objbias))
    sopts.objbias = opts.objbias;
end
if(isfield(opts,'display') && ~isempty(opts.display))
    sopts.display = dispLevel(opts.display);
end
sopts.optiver = optiver;

%Run SCIP
[x,fval,exitflag,stats] = scip(file,sopts);

%Assign Outputs
info.BBNodes = stats.BBnodes;
info.BBGap = stats.BBgap;
info.Time = toc(t);
info.Algorithm = 'SCIP: Spatial Branch and Bound using IPOPT and SoPlex [AMPL Interface]';

%Process Return Code
[info.Status,exitflag] = scipRetCode(exitflag);


function  print_level = dispLevel(lev)
%Return CLP compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 4;
    case 'final'
        print_level = 3;
end