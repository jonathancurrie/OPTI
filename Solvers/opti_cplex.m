function [x,fval,exitflag,info] = opti_cplex(H,f,A,rl,ru,lb,ub,xint,sos,qc,x0,opts)
%OPTI_CPLEX Solve a LP/MILP/QP/MIQP/QCQP/MIQCQP using CPLEX
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         qrl <= x'Q'x + l'x <= qru
%                                              lb <= x <= ub
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1} 
%
%   x = opti_cplex([],f,A,rl,ru,lb,ub,xint) solves a LP/MILP where f is the 
%   objective vector, A,rl,ru are the linear constraints, lb,ub are the
%   bounds and xint is a string of integer variables ('C', 'I', 'B').
%
%   x = opti_cplex(H,f,A,rl,ru,lb,ub,xint) solves a QP/MIQP where H is the
%   objective matrix, and the remainder of the arguments are as above.
%
%   x = opti_cplex(H,...,xint,sos) sos is a structure with fields type, 
%   index, weight for SOS.
%
%   x = opti_cplex(H,...,sos,qc) qc is structure with fields Q, l, qrl and
%   qru for quadratic constraints.
%
%   x = opti_cplex(H,...,qc,x0) supples an initial guess of the solution.
%
%   x = opti_cplex(H,f,...,x0,opts) uses opts to pass CPLEX parameters to 
%   the solver. For example to set the node limit:
%       opts.mip.limits.nodes.Cur = 1000;
%
%   [x,fval,exitflag,info] = opti_cplex(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 12, opts = []; end
if nargin < 11, x0 = []; end
if nargin < 10, qc = []; end
if nargin < 9, sos = []; end
if nargin < 8, xint = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to opti_cplex'); end

% Set Defaults (not returned from Cplex if infeasible)
x = []; fval = [];

%Check cplex exists
if(isempty(which('Cplex.p')))
    error('Cannot find Cplex the class! Please ensure the MATLAB Cplex interface is installed on your computer!');
end

%Build CPLEX Model
cp = Cplex('OPTI');
%Fill in common properties
cp.Model.obj = f;
if(isempty(lb))
    cp.Model.lb = -Inf(size(f));
else
    cp.Model.lb = lb;
end
if(isempty(ub))
    cp.Model.ub = Inf(size(f));
else
    cp.Model.ub = ub;
end
if(isempty(A))
    cp.Model.A = zeros(0,length(f));
    cp.Model.lhs = zeros(0,0);
    cp.Model.rhs = zeros(0,0);
else
    cp.Model.A = A;
    if(isempty(rl))
        cp.Model.lhs = -Inf(size(A,1),1);
    else
        cp.Model.lhs = rl;
    end
    if(isempty(ru))
        cp.Model.rhs = Inf(size(A,1),1);
    else
        cp.Model.rhs = ru;
    end
end
%Check and fill in remainder of optional args
if(~isempty(x0)),   cp.Start.x = x0; end
if(~isempty(H)),    cp.Model.Q = H; end
if(~isempty(xint)), cp.Model.ctype = upper(xint); end
if(~isempty(qc))
    %Check we have all fields
    if(~isstruct(qc) || ~isfield(qc,'Q') || ~isfield(qc,'l') || ~isfield(qc,'qrl') || ~isfield(qc,'qru'))
        error('Quadratic constraints must be a structure with fields ''Q'', ''l'', ''qrl'' and ''qru''');
    end
    if(iscell(qc.Q)) %multiple constraints
        if(length(qc.Q) ~= size(qc.l,2) || length(qc.qrl) ~= length(qc.Q) || length(qc.qrl) ~= length(qc.qru))
            error('Error with quadratic constraint dimensions');
        end
        for i = 1:length(qc.Q)
            addQC(cp,qc.Q{i},qc.l(:,i),qc.qrl(i),qc.qru(i));
        end
    else %single row constraint
        addQC(cp,qc.Q,qc.l,qc.qrl,qc.qru);
    end
end
if(~isempty(sos))
    %Check we have all fields
    if(~isstruct(sos) || ~isfield(sos,'type') || ~isfield(sos,'index') || ~isfield(sos,'weight'))
        error('SOS constraints must be a structure with fields ''type'', ''index'' and ''weight''');
    end
    %Ensure rows for Cplex
    if(iscell(sos.index) && size(sos.index,1) > size(sos.index,2))
        sos.index = sos.index';
    end
    if(iscell(sos.weight) && size(sos.weight,1) > size(sos.weight,2))
        sos.weight = sos.weight';
    end
    cp.addSOSs(sos.type,sos.index,sos.weight); 
end

%Get objective constant term
if(~isempty(opts) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    objbias = opts.objbias;
    opts = rmfield(opts,'objbias');
else
    objbias = 0;
end

%Assign User supplied options to Cplex
if(~isempty(opts))
    on = fieldnames(opts);
    if(length(on) == 1 && strcmpi('Param',on{1})) %remove Param if present
        opts = opts.Param;
    end
    %Bit of a hack to be able to set Cplex parameters, must be an easier way!
    fn = nestedfieldnames(opts); 
    %Get cell arrays of each field
    cn = regexp(fn,'\.','split');
    %For each field, assign to param struct
    param = cp.Param;
    for i = 1:length(cn)
        try
            param = setfield(param,cn{i}{:},getfield(opts,cn{i}{:}));
        catch ME
            error('Error setting a Cplex parameter [%s]. Ensure you are setting the .Cur field of a valid parameter!\n\nError: %s',fn{i},ME.message);
        end
    end
    %Reassign parameters
    cp.Param = param;
end
    
%OPTI Default Parameters
cp.Param.output.clonelog.Cur = 0;

%Solve problem
sol = cp.solve();

%Assign Outputs
if(isfield(sol,'x')), x = sol.x; end
if(isfield(sol,'objval')), fval = sol.objval; end
if(isfield(sol,'itcnt')), info.Iterations = sol.itcnt; end
if(isfield(sol,'mipitcnt')), info.Nodes = sol.mipitcnt; end
if(isfield(sol,'miprelgap')), info.RelGap = sol.miprelgap; end
info.Time = toc(t); ex = '';
if(isfield(sol,'method'))    
    switch(sol.method)
        case 0, ex = 'Automatic';
        case 1, ex = 'Primal Simplex';
        case 2, ex = 'Dual Simplex';
        case 3, ex = 'Network Simplex';
        case 4, ex = 'Barrier';
        case 5, ex = 'Sifting';
        case 6, ex = 'Concurrent Dual and Barrier';
        case 11, ex = 'Feasible Optimal';
        case 12, ex = 'MIP Solver';
        case 13, ex = 'Robust Algorithm';
    end
end
info.Algorithm = ['Cplex: ' ex];
if(isfield(sol,'statusstring'))
    info.Status = sol.statusstring;
    info.Status(1) = upper(info.Status(1));
else
    info.Status = 'Unknown';
end
if(isfield(sol,'status'))
    switch(sol.status)
        case {1,101,102,24}
            exitflag = 1; %ok
        case {25,10,11,12,105,107,111,112}
            exitflag = 0; %time/iterations/nodes
        case {3,5,103,115,106}
            exitflag = -1; %infeasible
        case {4,118,119}
            exitflag = -2; %unbounded
        case {13}
            exitflag = -5;
        otherwise
            exitflag = -6;
    end
else
    exitflag = -6;
end

%Add constant objective term
if(~isempty(fval))
    fval = fval + objbias;
end

%Assign Lambda
if(isfield(cp.Solution,'dual') && ~isempty(cp.Solution.dual))
    info.Lambda.lin = cp.Solution.dual; %not splitting up into ineq and eq as we are using row constraints
    rc = cp.Solution.reducedcost;
    zu = rc; zu(zu > 0) = 0;
    zl = rc; zl(zl < 0) = 0;
    info.Lambda.upper = zu;
    info.Lambda.lower = zl;
end


function addQC(cp,Q,l,rl,ru)
if(isinf(rl))
    cp.addQCs(l,Q,'L',ru); 
elseif(isinf(ru))
    cp.addQCs(l,Q,'G',rl);
elseif(rl == ru)
    error('Cplex does not support quadratic equalities!');
else
    cp.addQCs(l,Q,'L',ru);
    cp.addQCs(l,Q,'G',rl);
end
