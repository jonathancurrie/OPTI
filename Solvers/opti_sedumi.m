function [x,fval,exitflag,info] = opti_sedumi(At,b,c,K,opts)
%OPTI_SEDUMI Solve a SDP using SEDUMI
%
%   x = opti_sedumi(At,b,c,K) solves a LP or SDP in SeDuMi format.
%
%   x = opti_sedumi(At,b,c,K,opts) uses opts to pass optiset options to the
%   solver. 
%
%   [x,fval,exitflag,info] = opti_sedumi(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR SEDUMI (GPL)

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 5), opts = optiset; end
if(nargin < 4), K = []; end
if(nargin < 3), error('This function expects at least 3 arguments!'); end

%Set SEDUMI options
pars = opts.solverOpts;
pars.fid = min(dispLevel(opts.display),1);
pars.maxiter = opts.maxiter;
pars.eps = opts.tolrfun;

t = tic;
%Run Solver
[X,x,sinfo] = sedumi(At,b,c,K,pars);

%Assign primal fval
fval = -b'*x;

%Assign Outputs
info.Iterations = sinfo.iter;
info.Time = toc(t);
info.Algorithm = 'SeDuMi: Self-Dual-Minimization SDP Solver';

%Assign Exit Flag
if(sinfo.iter >= opts.maxiter)
    info.Status = 'Exceeded Maximum Iterations';
    exitflag = 0;
elseif(sinfo.numerr > 0)
    info.Status = 'Numerical Errors';
    exitflag = -2;
elseif(sinfo.pinf==1 || sinfo.dinf==1)
    info.Status = 'Infeasible';
    exitflag = -1;
elseif(sinfo.pinf==0 && sinfo.dinf==0)
    info.Status = 'Optimal';
    exitflag = 1;
end
info.DualObjective = c'*X;
info.X = X;

