function [x,fval,exitflag,info] = opti_lpsolve(f,A,b,Aeq,beq,lb,ub,int,sos,opts)
%OPTI_LPSOLVE Solve a LP or MILP using LP_SOLVE
%
%   [x,fval,exitflag,info] = opti_lpsolve(f,A,b,Aeq,beq,lb,ub,int,sos) solves 
%   the linear program min f'x where A,b are the inequality constraints, Aeq,
%   beq are the equality constraints, lb,ub are the bounds and int is a 
%   string of the integer variables in the form 'CIB' (Continuous / Integer 
%   Binary). sos is a structure with fields sostype, sosind, soswt for SOS.
%
%   THIS IS A WRAPPER FOR LP_SOLVE USING THE MEX INTERFACE
%   See supplied LGPL License

%   Copyright (C) 2011 Jonathan Currie (IPL)

% Modified Version of lp_solve.m supplied with LP_SOLVE distribution

t = tic;

% Handle missing arguments
if nargin < 10, opts = optiset; end
if nargin < 9, sos = []; end
if nargin < 8, int = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_lpsolve'); end

%Check for linear constraints & augment as required
if(~isempty(A) || ~isempty(Aeq)) 
    %All <= for A,b
    con_type = ones(1,size(A,1)); %1 <=
    %Augment A, Aeq
    if(~isempty(Aeq))
        A = [A;Aeq]; b = [b;beq];
        con_type = [con_type 3*ones(1,size(Aeq,1))]; %3 ==
    end
    %Warn re sparse
    if(~issparse(A))
        if(strcmpi(opts.warnings,'all'))
            warning('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
        end
        A = sparse(A);
    end
end

%Load LP_SOLVE Problem
n = length(f);
m = size(A,1);
lp = lp_solve('make_lp', m, n);
%Load Objective
lp_solve('set_obj_fn', lp, f);
lp_solve('set_minim', lp); % default is solving minimum lp.
%Load Linear Constraints
if(~isempty(A))
    lp_solve('set_mat', lp, A);
    lp_solve('set_rh_vec', lp, b);
    lp_solve('set_constr_type', lp, con_type);
end
%Load Bounds
if(~isempty(lb)), lp_solve('set_lowbo', lp, lb); end
if(~isempty(ub)), lp_solve('set_upbo', lp, ub); end
%Set Integer Constraints
ind = find(int ~= 'C');
ismip = 0;
if(any(ind))
    for i = 1:length(ind)
      if(int(ind(i)) == 'B')
          lp_solve('set_binary', lp, ind(i), 1);
      elseif(int(ind(i)) == 'I')
          lp_solve('set_int', lp, ind(i), 1);
      end
    end
    ismip = 1;
end
%Set SOS
if(~isempty(sos) && isfield(sos,'type') && ~isempty(sos.type))
    if(length(sos.type) > 1)
        [r,c] = size(sos.type);
        if(~iscell(sos.index) || ~iscell(sos.weight))
            error('When adding multiple SOS they must be supplied as cell arrays!');
        end
        if(r > c)
            tp = str2num(sos.type); %#ok<*ST2NM>
        else
            tp = str2num(sos.type');
        end
        for i = 1:length(t) %not sure if this is actually working?
            lp_solve('add_SOS', lp, ['set' num2str(i)], tp(i), i, sos.index{i}, sos.weight{i});
        end
    else
        lp_solve('add_SOS', lp, 'set1', str2num(sos.type), 1, sos.index, sos.weight);
    end
end
%Set Options
lp_solve('set_verbose', lp, dispLevel(opts.display)); 
lp_solve('set_timeout', lp, opts.maxtime); %converts internally to long?
lp_solve('set_bb_depthlimit', lp, opts.maxnodes); %not exactly number of nodes..
lp_solve('set_epsint', lp, opts.tolint);
%Hard Coded Options
lp_solve('set_scaling', lp, 7); %curtis-reid scaling
lp_solve('set_presolve', lp, 1+2+4); %presolve rows + cols + elim linear dependent rows
%Solve
exitflag = lp_solve('solve', lp);
%Get Outputs
[fval, x] = lp_solve('get_solution', lp);
iter = lp_solve('get_total_iter', lp);
msg = lp_solve('get_statustext', lp, exitflag);
if(~ismip)
    try
        lam = lp_solve('get_dual_solution',lp);
    catch %#ok<CTCH>
        lam = struct('ineqlin',[],'eqlin',[],'bounds',[]);
    end
end
%Delete Problem
lp_solve('delete_lp', lp);
%Check for deleted vars
if(length(x) ~= length(f))
    optiwarn('LP_SOLVE presolve has removed %d variable(s) from the optimization problem, thus the returned decision variable vector will be shorter than expected. OPTI does not currently insert the missing variables.',abs(length(x) ~= length(f)));
end

%Assign Outputs
info.Iterations = iter;
info.Time = toc(t);
info.Algorithm = 'LP_SOLVE: Simplex';
info.Status = msg;

switch(exitflag)
    case {0,1} %suboptimal as well...?
        exitflag = 1;
    case {2,3} %infeasible
        exitflag = -1;
    case 7 %timeout
        exitflag = 0;
    otherwise %unknown
        exitflag = -2;
end

%Add constant objective term
if(~isempty(fval) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    fval = fval + opts.objbias;
end

%Assign Lambda
if(~ismip)
    eq = con_type == 3;
    info.Lambda = struct('ineqlin',lam(~eq),'eqlin',lam(eq),'bounds',lam(length(b)+1:end));
else
    info.Lambda = [];
end

function  print_level = dispLevel(lev)
%Return LPSOLVE compatible display level
switch(lower(lev))
    case'off'
        print_level = 3;
    case 'iter'
        print_level = 6;
    case 'final'
        print_level = 4;
end