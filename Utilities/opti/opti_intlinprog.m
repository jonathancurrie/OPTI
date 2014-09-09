function [x,fval,exitflag,info,Opt] = opti_intlinprog(f,int,A,b,Aeq,beq,lb,ub,opts)
%OPTI_INTLINPROG Solve a MILP using an OPTI MILP Solver (Matlab-Like Overload)
%
%   [x,fval,exitflag,info] = opti_intlinprog(f,int,A,b,Aeq,beq,lb,ub) solves 
%   the linear program min f'x where A,b are the inequality constraints, 
%   Aeq,beq are the equality constraints, lb,ub are the decision variable 
%   bounds and int are the integer variable indices (double vector).
%
%   [x,fval,exitflag,info] = opti_intlinprog(f,...,ub,opts) allows the user 
%   to specify optiset options. This includes specifying a solver via the
%   'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_intlinprog(f,...) returns the internally built
%   OPTI object.

%   Copyright (C) 2014 Jonathan Currie (I2C2)


% Handle missing arguments
if nargin < 9, opts = optiset; end 
if nargin < 8, ub = []; end
if nargin < 7, lb = []; end
if nargin < 6, beq = []; end
if nargin < 5, Aeq = []; end
if nargin < 4, error('You must supply at least 4 arguments to opti_intlinprog'); end

%Build OPTI Object
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'int',int,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);
