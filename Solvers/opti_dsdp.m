function [x,fval,exitflag,info] = opti_dsdp(f,A,b,lb,ub,sdcone,x0,opts)
%OPTI_DSDP Solve a SDP using DSDP [Primal Form]
%
%   min f'*x      subject to:     A*x <= b
%    x                            lb <= x <= ub
%                                 X = F1*x1 + F2*x2 + ... + Fn*xn - F0
%                                 X >= 0 [positive semidefinite]
%                                 
%                                 
%
%   x = opti_dsdp(f,A,b,lb,ub) solves a LP where f is the objective 
%   vector, A,b are the linear inequality constraints, and lb,ub are the 
%   bounds.
%
%   x = opti_dsdp(f,...,ub,sdcone) uses the cell array sdcone to pass
%   semidefinite constraints to the solver. Each cell is a sparse matrix of
%   the form [F0(:) F1(:) F2(:)...] (or in alternative notation, [C A0 A1 A2..])
%   where each F matrix has been converted to a column vector and concatenated.
%
%   x = opti_dsdp(f,...,x0,opts) uses opts to pass optiset options to the
%   solver. 
%
%   [x,fval,exitflag,info] = opti_dsdp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR DSDP USING THE MEX INTERFACE
%   See supplied License

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

t = tic;

% Handle missing arguments
if nargin < 8, opts = optiset; end 
if nargin < 7, x0 = []; end 
if nargin < 6, sdcone = []; end 
if nargin < 5, ub = []; end
if nargin < 4, lb = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_dsdp'); end

warn = strcmpi(opts.warnings,'all');

%Check sparsity
if(~isempty(A) && ~issparse(A))
    if(warn), optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]'); end
    A = sparse(A);
end

%Convert SDP constraints into DSDP format
if(~isempty(sdcone))
    sdcone = DSDPSDCones(sdcone,length(f));
end

%Get objective constant term
if(~isempty(opts) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    objbias = opts.objbias;
    opts = rmfield(opts,'objbias');
else
    objbias = 0;
end

%Addin dsdp settings if specified
if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts))
    dopts = dsdpset(opts.solverOpts);    
else    
    dopts = [];
end
%Add OPTI Options
dopts.maxtime = opts.maxtime;
dopts.maxiter = opts.maxiter;
dopts.display = dispLevel(opts.display);
dopts.optiver = optiver;    

%MEX contains error checking
[x,fvals,exitflag,stats,X] = dsdp(-f, A, b, lb, ub, sdcone, x0, dopts);

%Assign X to matrices
if(~isempty(X))
    for i = 1:length(X)
        try
            X{i} = dmat(X{i});
        catch
        end
    end
end

%Invert fval
fval = -fvals.pval + objbias;

%Assign Outputs
info.Iterations = stats.iter;
info.Time = toc(t);
info.Algorithm = 'DSDP: Dual Scaling Interior Point SDP Solver';

switch(exitflag)
    case 1
        switch(stats.pdflag)
            case 0
                info.Status = 'Not sure whether Primal or Dual is feasible, check bounds';
                exitflag = -2;
            case 1
                info.Status = 'Both Primal and Dual are Feasible and Bounded';
            case 3
                info.Status = 'Dual is Unbounded, Primal Is Infeasible';
                exitflag = -1;
            case 4
                info.Status = 'Dual is Infeasible and Primal is Unbounded';
                exitflag = -1;
            otherwise
                info.Status = 'Unknown';
        end
    case -3
        info.Status = 'Exceeded Iterations';
        exitflag = 0;
    case -27
        info.Status = 'Exceeded Maximum Time';
        exitflag = 0;
    case -2
        info.Status = 'Failed: Numerical Difficulties with Short Step Lengths';
        exitflag = -3;
    case -6
        info.Status = 'Failed: Initial points imply S is not positive';
        exitflag = -3;
    case -8
        info.Status = 'Failed: Indefinite Schur Matrix';
        exitflag = -3;
    case -9
        info.Status = 'Failed: Numerical Error';
        exitflag = -3;
    case 5
        info.Status = 'Failed: Dual Objective big enough to stop';
        exitflag = -3;
    case 7
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise
        info.Status = [];
end

info.DualObjective = -fvals.dval + objbias;
info.X = X;

%Convert SDCone to DSDP Format (TRIU)
function sdout = DSDPSDCones(sdcone,n)
if(~iscell(sdcone))
    sdout = TriUCone(sdcone,n);    
else
    sdout = cell(length(sdcone),1);
    for i = 1:length(sdcone)
        sdout{i} = TriUCone(sdcone{i},n);
    end
end

%Guts of above
function dsdp_cone = TriUCone(cone,n)
m = sqrt(size(cone,1));
if((n+1) ~= size(cone,2)), error('Each F matrix must be supplied as a column vector (i.e. [F0(:) F1(:) ...]) such that there are n+1 columns'); end
if(floor(m) ~= m), error('Semidefinite constraint does not have a number of rows that can be converted into a square matrix'); end
ind = triu(ones(m))==1;
dsdp_cone = -cone(ind,:); %negate all to get to DSDP format

function [A]=dmat(V)
% DSDP5.0 
% Copyright (c) 2003 by
% S. Benson and Y. Ye
% Last modified: December 2003
[m,m1]=size(V);
if (m1~=1), error('Not a column vector.'); end;
n=floor(sqrt(2*m));
if (n*(n+1)/2 ~= m) 
     error('Impossible Dimension. Cannot convert to square matrix.'); 
end;
nzV=nnz(V);
if (issparse(V) || nzV < n*(n+1)/2), A=sparse(n,n); else A=zeros(n,n);  end;
for i=1:n, 
    k1=i*(i-1)/2+1;  k2=i*(i+1)/2;
    A(1:i,i)= V(k1:k2); 
end;
A=A+triu(A,1)';  



