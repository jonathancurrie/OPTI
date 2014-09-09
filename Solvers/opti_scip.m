function [x,fval,exitflag,info] = opti_scip(H,f,A,rl,ru,lb,ub,xint,sos,qc,opts)
%OPTI_SCIP Solve a LP/MILP/QP/MIQP/QCQP/MIQCQP using SCIP
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         qrl <= x'Qx + l'x <= qru
%                                              lb <= x <= ub
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1} 
%
%   x = opti_scip([],f,A,rl,ru,lb,ub,xint) solves a LP/MILP where f is the 
%   objective vector, A,rl,ru are the linear constraints, lb,ub are the
%   bounds and xint is a string of integer variables ('C', 'I', 'B').
%
%   x = opti_scip(H,f,A,rl,ru,lb,ub,xint) solves a QP/MIQP where H is the
%   objective matrix, and the remainder of the arguments are as above.
%
%   x = opti_scip(H,...,xint,sos) sos is a structure with fields type, 
%   index, and weight for SOS.
%
%   x = opti_scip(H,...,sos,qc) qc is structure with fields Q, l, and qrl 
%   and qru for quadratic constraints.
%
%   x = opti_scip(H,f,...,qc,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_scip(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR SCIP USING THE MEX INTERFACE
%   See supplied ZIB Academic License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 11, opts = optiset; end
if nargin < 10, qc = []; end
if nargin < 9, sos = []; end
if nargin < 8, xint = repmat('C',size(f)); end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_scip'); end

warn = strcmpi(opts.warnings,'all');

%Check sparsity
if(~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    A = sparse(A);
end
if(~isempty(H) && ~issparse(H))
    if(warn)
        optiwarn('opti:sparse','The H matrix should be sparse, correcting: [sparse(A)]');
    end
    H = sparse(H);
end
if(~isempty(qc))
    if(~isfield(qc,'Q') || ~isfield(qc,'l') || ~isfield(qc,'qrl') || ~isfield(qc,'qru'))
        error('The qc structure must contain fields Q, l, qrl and qru');
    end
    if(~isempty(qc.Q))
        err = 0;
        if(iscell(qc.Q))
            for i=1:length(qc.Q)
                if(~issparse(qc.Q{i}))
                    err = 1;
                    qc.Q{i} = sparse(qc.Q{i});
                end
            end                                       
        elseif(~issparse(qc.Q))
            err = 1;
            qc.Q = sparse(qc.Q);
        end
        if(err && warn)
            optiwarn('opti:sparse','The Q matrix should be sparse, correcting: [sparse(Q)]');
        end 
    end
end

%Remove H if all nz
if(~isempty(H) && nnz(H) == 0), H = []; end

%Setup Printing
opts.display = dispLevel(opts.display);

%MEX contains error checking
[x,fval,exitflag,stats] = scip(H, f, A, rl, ru, lb, ub, xint, sos, qc, [], opts);

%Assign Outputs
info.BBNodes = stats.BBnodes;
info.BBGap = stats.BBgap;
info.Time = toc(t);
if(~isempty(H) || ~isempty(qc))
    info.Algorithm = 'SCIP: Spatial Branch and Bound using IPOPT & SoPlex';
else
    info.Algorithm = 'SCIP: Spatial Branch and Bound using SoPlex';
end

switch(exitflag)
    case 10
        if(any(xint ~= 'C'))
            info.Status = 'Globally Integer Optimal';
        else
            info.Status = 'Globally Optimal';
        end
        exitflag = 1;
    case {2,3,4,5,6,7,8}
        info.Status = 'Exceeded Iterations / Time / Nodes';
        exitflag = 0;
    case 11
        info.Status = 'Infeasible';
        exitflag = -1;
    case {12,13}
        info.Status = 'Unbounded or Infeasible';
        exitflag = -2;
    case 1
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise
        info.Status = [];
end


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