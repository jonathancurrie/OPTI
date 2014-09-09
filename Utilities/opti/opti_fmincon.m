function [x,fval,exitflag,info,Opt] = opti_fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,opts)
%OPTI_FMINCON Solve a NLP using an OPTI NLP Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon) 
%   solves the constrained nonlinear optimization min f(x) where fun is 
%   the nonlinear function to be minimized [fun(x)], starting at x0. A,b
%   are linear inequality constraints, Aeq,beq are linear equality
%   constraints, lb,ub are decision variable bounds and nonlcon are the
%   nonlinear constraints in Matlab form.
%
%   [x,fval,exitflag,info] = opti_fmincon(fun,...,opts) allows the 
%   user to specify optiset options. This includes specifying a solver via 
%   the 'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_fmincon(fun,...) returns the internally 
%   built OPTI object.

%   Copyright (C) 2011 Jonathan Currie (I2C2)


%Handle missing arguments
if nargin < 10, opts = optiset; end 
if nargin < 9, nonlcon = []; end
if nargin < 8, ub = []; end
if nargin < 7, lb = []; end
if nargin < 6, beq = []; end
if nargin < 5, Aeq = []; end
if nargin < 4, b = []; end
if nargin < 3, A = []; end
if nargin < 2, error('You must supply at least 2 arguments to opti_fmincon'); end

%Sort out Fun + Grad
[f,g] = detGrad(fun,x0);

%Sort out NLCON
[nlcon,nlrhs,nle] = detNlcon(nonlcon,x0);

%Build OPTI Object
Opt = opti('fun',f,'grad',g,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'x0',x0,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);