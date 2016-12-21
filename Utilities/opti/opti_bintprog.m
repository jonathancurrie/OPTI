function [x,fval,exitflag,info,Opt] = opti_bintprog(f,A,b,Aeq,beq,x0,opts)
%OPTI_BINTPROG Solve a BIP using an OPTI BIP Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_bintprog(f,A,b,Aeq,beq,x0) solves the 
%   linear program min f'x where A,b are the inequality constraints, 
%   Aeq,beq are the equality constraints and x0 is the initial guess. All
%   decision variables are binary only.
%
%   [x,fval,exitflag,info] = opti_bintprog(f,...,x0,opts) allows the user 
%   to specify optiset options. This includes specifying a solver via the
%   'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_bintprog(f,...) returns the internally built
%   OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 7, opts = optiset; end 
if nargin < 6, x0 = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_bintprog'); end

%Build OPTI Object
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'x0',x0,'int',repmat('B',size(f)),'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);
