function [x,fval,exitflag,info,lambda,Opt] = opti_linprog(f,A,b,Aeq,beq,lb,ub,x0,opts)
%OPTI_LINPROG Solve a LP using an OPTI LP Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_linprog(f,A,b,Aeq,beq,lb,ub,x0) solves the 
%   linear program min f'x where A,b are the inequality constraints, 
%   Aeq,beq are the equality constraints and lb,ub are the bounds.
%
%   [x,fval,exitflag,info] = opti_linprog(f,...,ub,opts) allows the user to
%   specify optiset options. This includes specifying a solver via the
%   'solver' field of optiset.
%
%   [x,...,info,lambda] = opti_linprog(f,...) returns a structure of the 
%   Lagrange multipliers (dual solution).
%
%   [x,...,lambda,Opt] = opti_linprog(f,...) returns the internally built
%   OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 9, opts = optiset; end 
if nargin < 8, x0 = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_linprog'); end

%Build OPTI Object
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'x0',x0,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);

%Extra Dual Solution
if(isfield(info,'Lambda'))
    lambda = info.Lambda;
elseif(isfield(info,'X'))
    if(iscell(info.X))
        lambda.dual = full(info.X{1});
    else
        lambda.dual = full(info.X);
    end
else
    lambda = [];
end