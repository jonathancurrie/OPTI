function [x,fval,exitflag,info,Opt] = opti_mintprog(f,A,b,Aeq,beq,lb,ub,int,sos,opts)
%OPTI_MINTPROG Solve a MILP using an OPTI MILP Solver (Matlab-Like Overload)
%
%   [x,fval,exitflag,info] = opti_mintprog(f,A,b,Aeq,beq,lb,ub,int) solves 
%   the linear program min f'x where A,b are the inequality constraints, 
%   Aeq,beq are the equality constraints, lb,ub are the decision variable 
%   bounds and int is a string of integer variables ('C', 'I', 'B').
%
%   [x,fval,exitflag,info] = opti_minprog(f,...,int,sos) allows Special
%   Ordered Sets (SOS) to be specified via the sos structure. Fields must
%   include 'type' ('1' or '2'), 'index' (vector) and 'weight' (vector).
%   Multiple SOS can be supplied as cell arrays.
%
%   [x,fval,exitflag,info] = opti_mintprog(f,...,sos,opts) allows the user 
%   to specify optiset options. This includes specifying a solver via the
%   'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_mintprog(f,...) returns the internally built
%   OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 10, opts = optiset; end 
if nargin < 9, sos = []; end
if nargin < 8, int = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, error('You must supply at least 3 arguments to opti_mintprog'); end

%Check SOS
if(~isempty(sos))
    if(~isfield(sos,'type')), error('The SOS structure must have the field type!'); end
    if(~isfield(sos,'index')), error('The SOS structure must have the field index!'); end
    if(~isfield(sos,'weight')), error('The SOS structure must have the field weight!'); end
end

%Build OPTI Object
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'int',int,'sos',sos,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);
