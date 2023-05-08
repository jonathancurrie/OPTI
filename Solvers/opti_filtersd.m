function [x,fval,exitflag,info] = opti_filtersd(fun,grad,x0,lb,ub,nlcon,nljac,nljacstr,cl,cu,opts)
%OPTI_FILTERSD Solve a NLP using FILTERSD 
%
%   min f(x)       subject to:    cl <= nlcon(x) <= cu
%    x                            lb <= x <= ub
%
%   x = opti_filtersd(fun,grad,x0) solves a NLP where fun is the objective 
%   function, grad is the gradient of the objective and x0 is a starting 
%   guess.
%
%   x = opti_filtersd(fun,grad,x0,lb,ub) solves subject to decision 
%   variable bounds lb <= x <= ub. Infinite bounds are OK.
%
%   x = opti_filtersd(fun,x0,lb,ub,nlcon,nljac,[],cl,cu) solves subject to 
%   the nonlinear inequality constraints cl <= nlcon(x) <= cu, and nljac
%   returns a DENSE Jacobian of the constraints.
%
%   x = opti_filtersd(fun,x0,lb,ub,nlcon,nljac,nljacstr,cl,cu) solves as
%   above, except nljac returns a SPARSE Jacobian of the constraints, and
%   nljacstr supplies a sparse matrix of ALL possible non-zero locations
%   within the Jacobian.
%
%   x = opti_filtersd(fun,...,cu,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_filtersd(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR FILTERSD
%   See supplied Eclipse Public License

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

if(nargin < 11), opts = optiset; else opts = optiset(opts); end
if(nargin < 10), cu = []; end
if(nargin < 9),  cl = []; end
if(nargin < 8),  nljacstr = []; end
if(nargin < 7),  nljac = []; end
if(nargin < 6),  nlcon = []; end
if(nargin < 5),  ub = []; end
if(nargin < 4),  lb = []; end
if(nargin < 3),  error('FILTERSD requires at least 3 arguments'); end

%Default is not sparse
sp = false;

%Determine we have a sparse problem
if(~isempty(nljacstr))
    sp = true; %must be sparse if we passed an arg here
    if(~isa(nljacstr,'function_handle'))
        error('Jacobian Structure must be a function handle');
    end
end

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('FILTERSD requires an initial guess, x0!');
end

%Setup display level
opts.display = dispLevel(opts.display);
opts.optiver = optiver;

t = tic;
if(sp)
    % Run FILTERSD [Sparse]
    jacstr = nljacstr(); %filtersd just requires the matrix
    [x, fval, exitflag, stats, lambda] = filtersdsp(fun,grad,x0,lb,ub,nlcon,nljac,jacstr,cl,cu,opts);
else
    % Run FILTERSD [Dense]
    [x, fval, exitflag, stats, lambda] = filtersd(fun,grad,x0,lb,ub,nlcon,nljac,cl,cu,opts);
end

%Collect Results
info.Iterations = stats.niter;
info.FuncEvals = stats.nfval;
info.GradEvals = stats.ngval;
info.Time = toc(t);
if(sp)
    info.Algorithm = 'FILTERSD: Nonlinear Optimization using the Filter Search & Trust Region [Sparse]';
else
    info.Algorithm = 'FILTERSD: Nonlinear Optimization using the Filter Search & Trust Region [Dense]';
end

switch(exitflag)
    case 0
        info.Status = 'Locally Optimal';
        exitflag = 1;
    case {5,101,102}
        info.Status = 'Exceeded Iterations / Function Evaluations / Time';
        exitflag = 0;
    case {1,2,4,7,8}
        info.Status = 'Initialization Error';
        exitflag = -2;
    case {3,6}
        info.Status = 'Infeasible / Inconsistent Constraints';
        exitflag = -1;
    case {9,10}
        info.Status = 'LCP Solver Error';
        exitflag = -3;
    case 105
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise        
        info.Status = 'FILTERSD Error';
end

info.Lambda = lambda;
