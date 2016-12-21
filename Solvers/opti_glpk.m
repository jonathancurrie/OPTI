function [x,fval,exitflag,info] = opti_glpk(f,A,b,Aeq,beq,lb,ub,int,opts)
%OPTI_GLPK Solve a LP using the GLPK
%
%   [x,fval,exitflag,info] = opti_glpk(f,A,b,Aeq,beq,lb,ub,int) solves the 
%   linear program min f'x where A,b are the inequality constraints, Aeq,beq 
%   are the equality constraints and lb,ub are the bounds. int is a string of 
%   the integer variables in the form 'CIB' (Continuous / Integer / Binary)
%
%   THIS IS A WRAPPER FOR GLPK USING GLPKMEX
%   See supplied GNU License

%   Copyright (C) 2011 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 9, opts = optiset; end 
if nargin < 8, int = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3
    disp('GLPK Matlab interface. Version: 2.10');
    disp('(C) 2001-2007, Nicolo'' Giorgetti.');
    disp('Maintained by Niels Klitgord');
    disp(' ')
    disp('MATLAB front end Modified for OPTI Toolbox');
    return;
end

%All <= for A,b
ctype = char('U' * ones(1,size(A,1)));
%Augment Aeq
if(~isempty(Aeq))
    A = [A;Aeq]; b = [b;beq];
    ctype = [ctype char('S' * ones(1,size(Aeq,1)))];
end

%Ensure f is a column
f = f(:); nx = length(f);
%Make sure we have A
if (isempty(A))
    error('A cannot be an empty matrix');
end
[nc, nxa] = size(A);
if (~isreal(A) || nxa ~= nx)
    error('A must be a real valued %d by %d matrix', nc, nx);
end
if(~issparse(A))
    if(strcmpi(opts.warnings,'all'))
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    A = sparse(A);
end
%Make sure we have b
if (isempty(b))
    error('B cannot be an empty vector');
end
if (~isreal(b) || length(b) ~= nc)
    error('B must be a real valued %d by 1 vector', nc);
end
%Lower Bounds
if (isempty(lb))
    lb = -Inf(nx, 1);
elseif (~isreal(lb) || all(size(lb) > 1) || length(lb) ~= nx)
    error('LB must be a real valued %d by 1 column vector', nx);
end
%Upper Bounds
if (isempty(ub))
    ub = Inf(nx, 1);
elseif (~isreal(ub) || all(size(ub) > 1) || length(ub) ~= nx)
    error('UB must be a real valued %d by 1 column vector', nx);
end
%Binary / Integer
if isempty(int)
    int = repmat('C', nx, 1);
elseif (~ischar(int) || all(size(int) > 1) || length (int) ~= nx)
    error('INT must be a char valued vector of length %d', nx);
else
   for i=1:length(int)
      switch(int(i))
         case {'c','C'}, % do nothing
         case {'i','I'}, % do nothing
         case {'b','B'}, % do nothing
         otherwise
            error('VARTYPE must contain only C, I or B');
      end
   end
end

%Setup Printing
param.msglev = dispLevel(opts.display);
param.outfrq = 10;
%Setup Other Options
param.itlim = opts.maxiter;
param.tmlim = opts.maxtime*60; %I think this is in ms?
param.tolobj = opts.tolrfun;
param.tolint = opts.tolint;

%Solve Problem
[x,fval,exitflag,extra] = glpk(f, A, b, lb, ub, ctype, int, 1, param);

%Assign Outputs
info.Nodes = [];
info.Time = toc(t);
info.Algorithm = 'GLPK: Revised Simplex';

switch(exitflag)
    case 5
        stat = 'optimal';
        exitflag = 1;
    case {107,108}
        stat = 'iterations exceeded';
        exitflag = 0;
    case 6
        stat = 'unbounded';
        exitflag = -1;
    otherwise
        stat = 'infeasible /  error';
        exitflag = -1;
end

%Add constant objective term
if(~isempty(fval) && isfield(opts,'objbias') && ~isempty(opts.objbias))
    fval = fval + opts.objbias;
end

info.Status = stat;  
%Assign Lambda
if(isempty(int) || all(int == 'C'))
    eq = ctype == 'S';
    info.Lambda = struct('ineqlin',extra.lambda(~eq),'eqlin',extra.lambda(eq),'bounds',extra.redcosts);
else
    info.Lambda = [];
end

function  print_level = dispLevel(lev)
%Return GLPK compatible display level
switch(lower(lev))
    case'off'
        print_level = 1;
    case 'iter'
        print_level = 3;
    case 'final'
        print_level = 3;
end