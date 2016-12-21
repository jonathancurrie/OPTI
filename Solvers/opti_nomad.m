function [x,fval,exitflag,info] = opti_nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts)
%OPTI_NOMAD Solve a Global MINLP / NLP using NOMAD (Nonlinear Optimization using the MADS Algorithm)
%
%   min f(x)       subject to:    nlcon(x) <= nlrhs
%    x                            lb <= x <= ub
%                                 xi in Z
%                                 xj in {0,1} 
%
%   x = opti_nomad(fun,x0) solves a Global NLP where fun is the 
%   objective function and x0 is a starting guess.
%
%   x = opti_nomad(fun,x0,lb,ub) solves subject to decision variable bounds
%   lb <= x <= ub. Infinite bounds will be ignored by the solver.
%
%   x = opti_nomad(fun,x0,lb,ub,nlcon,nlrhs) solves subject to the
%   nonlinear inequality constraints nlcon(x) <= nlrhs.
%
%   x = opti_nomad(fun,...,nlrhs,xtype) solves subject to integer and
%   binary constraints specified in xtype. xtype is a string of characters
%   ('C', 'I', 'B') for continuous, integer, or binary constraints.
%
%   x = opti_nomad(fun,...,xtype,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_nomad(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR NOMAD
%   See supplied Lesser GNU Public License

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(nargin < 8), opts = optiset; else opts = optiset(opts); end
if(nargin < 7), xtype = []; end
if(nargin < 6), nlrhs = []; end
if(nargin < 5), nlcon = []; end
if(nargin < 4), ub = []; end
if(nargin < 3), lb = []; end
if(nargin < 2), error('NOMAD requires at least 2 arguments'); end

%Addin nomad settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    nopts = nomadset(opts.solverOpts);    
else    
    nopts = [];
end
%Override with optiset settings (note NOMAD does not accept unknown (optiset) settings)
nopts = nomadset(nopts,'max_iterations',opts.maxiter,'max_bb_eval',opts.maxfeval,...
                'max_time',opts.maxtime,'display_degree',dispLevel(opts.display),...
                'min_mesh_size',num2str(opts.tolrfun),'iterfun',opts.iterfun);

%If we have any inf bounds, set initial mesh size (C.Tribes idea)
if(any(isinf(lb)) || any(isinf(ub)) && isempty(nopts.initial_mesh_size))
    nopts = nomadset(nopts,'initial_mesh_size','10'); %better heuristic on size here?
end
%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('NOMAD requires an initial guess, x0!');
end
            
t = tic;
% Run NOMAD
[x, fval, exitflag, iter, nfval] = nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,nopts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = nfval;
info.Time = toc(t);
info.Algorithm = 'NOMAD: Nonlinear Optimization using the MADS Algorithm';

switch(exitflag)
    case 1
        info.Status = 'Converged / Target Reached';
    case 0
        info.Status = 'Exceeded Iterations / Function Evaluations / Time';
    case -1
        info.Status = 'Infeasible / Mesh Limit Reached';
    case -2
        info.Status = 'Initialization Error';
    case -3
        info.Status = 'NOMAD Error';
    case -5
        info.Status = 'User Exited';
    otherwise        
        info.Status = 'NOMAD Error';
end
