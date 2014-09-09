function [x,fval,exitflag,info] = opti_csdp(f,A,b,lb,ub,sdcone,x0,opts)
%OPTI_CSDP Solve a SDP using CSDP [Primal Form]
%
%   min f'*x      subject to:     A*x <= b
%    x                            lb <= x <= ub
%                                 X = F1*x1 + F2*x2 + ... + Fn*xn - F0
%                                 X >= 0 [positive semidefinite]
%                                 
%                                 
%
%   x = opti_csdp(f,A,b,lb,ub) solves a LP where f is the objective 
%   vector, A,b are the linear inequality constraints, and lb,ub are the 
%   bounds.
%
%   x = opti_csdp(f,...,ub,sdcone) uses the cell array sdcone to pass
%   semidefinite constraints to the solver. Each cell is a sparse matrix of
%   the form [F0(:) F1(:) F2(:)...] (or in alternative notation, [C A0 A1 A2..])
%   where each F matrix has been converted to a column vector and concatenated.
%
%   x = opti_csdp(f,...,x0,opts) uses opts to pass optiset options to the
%   solver. 
%
%   [x,fval,exitflag,info] = opti_csdp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR CSDP USING THE MEX INTERFACE
%   See supplied Eclipse Public License

%   Copyright (C) 2013 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end 
if nargin < 7, x0 = []; end 
if nargin < 6, sdcone = []; end 
if nargin < 5, ub = []; end
if nargin < 4, lb = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_csdp'); end

warn = strcmpi(opts.warnings,'all');

%Check sparsity
if(~isempty(A) && ~issparse(A))
    if(warn), optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]'); end
    A = sparse(A);
end
if(~isempty(sdcone))
    ok = 1;
    if(iscell(sdcone))        
        for i = 1:length(sdcone)
            if(~issparse(sdcone{i}))
                sdcone{i} = sparse(sdcone{i});
                ok = 0;
            end
        end
    elseif(~issparse(sdcone))        
        sdcone = sparse(sdcone);
        ok = 0;
    end
    if(~ok && warn)
        optiwarn('opti:sparse','The Semidefinite Constraints should be sparse, correcting: [sparse(sdcone)]');
    end
end

%Addin csdp settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    copts = csdpset(opts.solverOpts);    
else    
    copts = [];
end
%Add OPTI Options
copts.maxtime = opts.maxtime;
copts.maxiter = opts.maxiter;
copts.display = dispLevel(opts.display);
    
%MEX contains error checking
[x,fvals,exitflag,stats,X] = csdp(f, A, b, lb, ub, sdcone, x0, copts);

%Assign primal fval
fval = fvals.pval;

%Assign Outputs
info.Iterations = stats.iter;
info.Time = toc(t);
info.Algorithm = 'CSDP: Predictor-Corrector Primal-Dual SDP Solver';

switch(exitflag)
    case 0
        info.Status = 'CSDP Converged';
        exitflag = 1;
    case 1
        info.Status = 'Primal Infeasible';
        exitflag = -1;
    case 2
        info.Status = 'Dual Infeasible';
        exitflag = -1;
    case 3
        info.Status = 'Partial Success (Full Accuracy Not Achieved)';
        exitflag = 3;
    case 4
        info.Status = 'Exceeded Maximum Iterations';
        exitflag = 0;
    case 5
        info.Status = 'Stuck at Edge of Primal Feasibility';
        exitflag = -2;
    case 6
        info.Status = 'Stuck at Edge of Dual Infeasibility';
        exitflag = -2;
    case 7
        info.Status = 'Lack of Progress';
        exitflag = -3;
    case 8
        info.Status = 'X, Z or O was singular';
        exitflag = -3;
    case 9
        info.Status = 'Detected NaN or Inf';
        exitflag = -3;
    case 10
        info.Status = 'General easy_sdp() failure';
        exitflag = -3;
    case 11
        info.Status = 'Failed check on C (check symmetry)';
        exitflag = -3;
    case 12
        info.Status = 'Failed constraint check';
        exitflag = -3;
    case -27
        info.Status = 'Exceeded Maximum Time';
        exitflag = 0;
    case -50
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise
        info.Status = [];
end

info.DualObjective = fvals.dval;
info.X = X;



