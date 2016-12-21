function [x,fval,exitflag,info,lambda,Opt] = opti_quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opts)
%OPTI_QUADPROG Solve a QP using an OPTI QP Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_quadprog(H,f,A,b,Aeq,beq,lb,ub,x0) solves 
%   the quadratic program min 1/2x'Hx + f'x where A,b are the inequality 
%   constraints, Aeq,beq are the equality constraints and lb,ub are the 
%   bounds.
%
%   [x,fval,exitflag,info] = opti_quadprog(H,...,ub,opts) allows the user 
%   to specify optiset options. This includes specifying a solver via the
%   'solver' field of optiset.
%
%   [x,...,info,lambda] = opti_quadprog(H,...) returns a structure of the 
%   Lagrange multipliers (dual solution).
%
%   [x,...,lambda,Opt] = opti_quadprog(H,...) returns the internally built
%   OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 10, opts = optiset; end 
if nargin < 9, x0 = []; end
if nargin < 8, ub = []; end
if nargin < 7, lb = []; end
if nargin < 6, beq = []; end
if nargin < 5, Aeq = []; end
if nargin < 4, error('You must supply at least 4 arguments to opti_quadprog'); end

%Build OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'x0',x0,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);

%Extra Dual Solution
if(isfield(info,'Lambda'))
    lambda = info.Lambda;
else
    lambda = [];
end
