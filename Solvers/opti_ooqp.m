function [x,fval,exitflag,info] = opti_ooqp(H,f,A,rl,ru,Aeq,beq,lb,ub,opts)
%OPTI_OOQP Solve a LP or QP using OOQP
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru (ineq only)
%    x                                         Aeq*x = beq
%                                              lb <= x <= ub
%
%   x = opti_ooqp(H,f,A,rl,ru,Aeq,beq,lb,ub) solves a QP where H and f are 
%   the objective matrix and vector respectively, A,rl,ru are the linear 
%   inequality constraints, Aeq,beq are the linear equality constraints and 
%   lb,ub are the bounds. Set H to [] for solving an LP.
%
%   x = opti_ooqp(H,...,ub,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_ooqp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR OOQP USING OOQP MEX INTERFACE
%   See supplied License

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_ooqp'); end

warn = strcmpi(opts.warnings,'all');

%Add in ooqpset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = ooqpset(opts.solverOpts);
    %Add in options from optiset (not all currently used)    
    popts.maxiter = opts.maxiter;
    popts.maxtime = opts.maxtime;
    popts.tolrfun = opts.tolrfun;
    %Add objective bias
    if(isfield(opts,'objbias')), popts.objbias = opts.objbias; end
else
    popts = opts;
end
%Setup display level
popts.display = dispLevel(opts.display);
popts.optiver = optiver;

%Ensure Sparsity
if(~isempty(A) && ~issparse(A))
    if(warn), optiwarn('OPTI:NotSparseA','The A matrix should be sparse, correcting: [sparse(A)]'); end
    A = sparse(A);
end
if(~isempty(H) && ~issparse(H))
    if(warn), optiwarn('OPTI:NotSparseH','The H matrix should be sparse, correcting: [sparse(H)]'); end
    H = sparse(H);
end
%Ensure Sym TRIU
if(any(any(triu(H,1) ~= 0)))
    if(warn), optiwarn('OPTI:NotTrilH','The H matrix should be Symmetric Lower Triangular, correcting: [tril(H)]'); end
    H = tril(H);
end

%Call Solver
[x,fval,stat,iter,lam] = ooqp( H, f, A, rl, ru, Aeq, beq, lb, ub, popts);

%Assign Outputs
info.Iterations = iter;
info.Time = toc(t);
if(~isfield(popts,'algorithm') || isempty(popts.algorithm))
    info.Algorithm = 'OOQP: Gondzio Predictor-Corrector with MA57';
else
    switch(popts.algorithm)
        case 0
            switch(popts.linear_solver)
                case 0
                    info.Algorithm = 'OOQP: Mehrotra Predictor-Corrector with PARDISO';
                case 1
                    info.Algorithm = 'OOQP: Mehrotra Predictor-Corrector with MA57';
                case 2
                    info.Algorithm = 'OOQP: Mehrotra Predictor-Corrector with MA27';
            end            
        case 1
            switch(popts.linear_solver)
                case 0
                    info.Algorithm = 'OOQP: Gondzio Predictor-Corrector with PARDISO';
                case 1
                    info.Algorithm = 'OOQP: Gondzio Predictor-Corrector with MA57';
                case 2
                    info.Algorithm = 'OOQP: Gondzio Predictor-Corrector with MA27';
            end
    end
end
            
switch(stat)
    case 0
        exitflag = 1;
        msg = 'Optimal';
    case 2
        exitflag = 0;
        msg = 'Exceeded Iterations';
    case 3 
        exitflag = -1;
        msg = 'Infeasible';
    case -1
        exitflag = 0;
        msg = 'Exceeded Maximum Time';
    case -2
        exitflag = -5;
        msg = 'User Exited';
    otherwise 
        exitflag = -2;
        msg = 'OOQP Error';
end

info.Status = msg;

%Assign Lambda
info.Lambda = struct('ineqlin',lam.pi,'eqlin',lam.y,'upper',lam.phi,'lower',lam.gamma);
    
    
    
    
    
