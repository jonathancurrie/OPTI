function [x,fval,exitflag,info] = opti_cbc(H,f,A,rl,ru,lb,ub,xtype,sos,x0,opts)
%OPTI_CBC Solve a MILP or MIQP using CBC
%
%   min f'*x      subject to:     rl <= A*x <= ru
%    x                            lb <= x <= ub
%                                 for i = 1..n: xi in Z
%                                 for j = 1..m: xj in {0,1} 
%
%   x = opti_cbc([],f,A,rl,ru,lb,ub,xtype) solves a MILP where f is the 
%   objective vector, A,rl,ru are the linear constraints, lb,ub are the
%   bounds and xtype is a string of integer variables ('C', 'I', 'B')
%
%   x = opti_cbc(H,f,A,rl,ru,lb,ub,xtype) solves a MIQP where H and f are 
%   the objective matrix and vector respectively, A,rl,ru are the linear 
%   constraints, lb,ub are the bounds and xtype is a string of integer 
%   variables ('C', 'I', 'B') [NOT CURRENTLY SUPPORTED]
%
%   x = opti_cbc([H,...,xtype,sos) sos is a structure with fields type, 
%   index, weight for SOS.
%
%   x = opti_cbc(H,...,sos,x0) uses x0 to warm start the solver.
%
%   x = opti_cbc(H,...,x0,opts) uses opts to pass cbcset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_cbc(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR CBC USING THE MEX INTERFACE
%   See supplied Eclipse Public License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 11, opts = optiset; end 
if nargin < 10, x0 = []; end
if nargin < 9, sos = []; end
if nargin < 8, xtype = repmat('C',size(f)); end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_cbc'); end

warn = strcmpi(opts.warnings,'all');

%Add in cbcset settings
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    popts = cbcset(opts.solverOpts);    
else
    popts = cbcset;
end
%Add in options from optiset   
popts.maxnodes = opts.maxnodes;
popts.maxtime = opts.maxtime;
popts.primalTol = opts.tolrfun;
popts.dualTol = opts.tolrfun;
popts.intTol = opts.tolint;
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
%Check SOS
if(~isempty(sos))
    if(~isstruct(sos) || ~isfield(sos,'type') || ~isfield(sos,'index') || ~isfield(sos,'weight'))
        error('SOS constraints must be a structure with fields ''type'', ''index'' and ''weight''');
    end
    if(length(sos.type) == 1)
        sos.index = {sos.index};
        sos.weight = {sos.weight};
    end
end

%MEX contains error checking
[x,fval,status,iter,cobj] = cbc(H, f, A, rl, ru, lb, ub, xtype, sos, x0, popts);

%Assign Outputs
info.Nodes = iter;
info.AbsGap = abs(cobj-fval);
info.RelGap = abs(cobj-fval)/(1e-1 + abs(fval));
info.Time = toc(t);
info.Algorithm = 'CBC: Branch and Cut using CLP';

switch(status)
    case 0
        info.Status = 'Integer Optimal';
        exitflag = 1;
    case 1
        info.Status = 'Linear Relaxation Infeasible';
        exitflag = -1;
    case 2
        info.Status = 'Gap Reached';
        exitflag = 1;
    case 3        
        info.Status = 'Maximum Nodes Reached';
        exitflag = 0;
    case 4
        info.Status = 'Maximum Time Reached';
        exitflag = 0;
    case 5
        info.Status = 'User Exited';
        exitflag = -5;
    case 6
        info.Status = 'Number of Solutions Reached';
        exitflag = 0;
    case 7
        info.Status = 'Linear Relaxation Unbounded';
        exitflag = -2;
	case 8
        info.Status = 'Proven Infeasible';
        exitflag = -1;
    otherwise
        info.Status = 'Unknown Termination';
        exitflag = -3;
end