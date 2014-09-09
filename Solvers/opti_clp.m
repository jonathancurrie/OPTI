function [x,fval,exitflag,info] = opti_clp(H,f,A,rl,ru,lb,ub,opts)
%OPTI_CLP Solve a LP or QP using CLP
%
%   min f'*x                   subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                 
%   OR
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%
%
%   x = opti_clp([],f,A,rl,ru,lb,ub) solves a LP where f is the objective 
%   vector, A,rl,ru are the linear constraints and lb,ub are the bounds.
%
%   x = opti_clp(H,f,A,rl,ru,lb,ub) solves a QP where H and f are 
%   the objective matrix and vector respectively, A,rl,ru are the linear 
%   constraints, and lb,ub are the bounds. 
%
%   x = opti_clp(H,...,ub,opts) uses opts to pass optiset options to the
%   solver. 
%
%   [x,fval,exitflag,info] = opti_clp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR CLP USING THE MEX INTERFACE
%   See supplied Eclipse Public License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end 
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_clp'); end

warn = strcmpi(opts.warnings,'all');

%Add in clpset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = clpset(opts.solverOpts);    
else
    popts = clpset;
end
%Add in options from optiset   
popts.maxiter = opts.maxiter;
popts.maxtime = opts.maxtime;
popts.primalTol = opts.tolrfun;
popts.dualTol = opts.tolrfun;
%Add objective bias
if(isfield(opts,'objbias')), popts.objbias = opts.objbias; end
%Setup display level
popts.display = dispLevel(opts.display,1,1);

%Check sparsity
if(~issparse(A))
    if(warn), optiwarn('OPTI:NotSparseA','The A matrix should be sparse, correcting: [sparse(A)]'); end
    A = sparse(A);
end
if(~issparse(H))
    if(~isempty(H) && warn), optiwarn('OPTI:NotSparseH','The H matrix should be sparse, correcting: [sparse(H)]'); end
    H = sparse(H);
end
%Check Sym Tril
if(any(any(triu(H,1) ~= 0)))
    if(warn), optiwarn('OPTI:NotTrilH','The H matrix should be Symmetric TRIL, correcting: [tril(H)]'); end
    H = tril(H);
end

%MEX contains error checking
[x,fval,exitflag,iter,lam] = clp(H, f, A, rl, ru, lb, ub, popts);

%Assign Outputs
info.Iterations = iter;
info.Time = toc(t);
switch(lower(popts.algorithm))
    case 'dualsimplex'
        info.Algorithm = 'CLP: Dual Simplex';
    case 'primalsimplex'
        info.Algorithm = 'CLP: Primal Simplex';
    case 'primalsimplexorsprint'
        info.Algorithm = 'CLP: Primal Simplex or Sprint';
    case 'barrier'
        info.Algorithm = 'CLP: Barrier';
    case 'barriernocross'
        info.Algorithm = 'CLP: Barrier (No Crossover)';
    case 'automatic'
        info.Algorithm = 'CLP: Automatically Chosen Solver';
end

switch(exitflag)
    case 0
        info.Status = 'Proven Optimal';
        exitflag = 1;
    case 1
        info.Status = 'Proven Primal Infeasible';
        exitflag = -1;
    case 2
        info.Status = 'Proven Dual Infeasible';
        exitflag = -1;
    case 3
        info.Status = 'Exceeded Maximum Iterations or Time';
        exitflag = 0;
    case 4
        info.Status = 'Clp Error';
        exitflag = -3;
    case 5
        info.Status = 'User Exit';
        exitflag = -5;
    otherwise
        info.Status = 'Unknown Exit';
        exitflag = -4;
end

%Assign Lambda
eq = rl == ru;
info.Lambda = struct('ineqlin',lam.dual_row(~eq),'eqlin',lam.dual_row(eq),'bounds',lam.dual_col);