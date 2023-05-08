function [x,fval,exitflag,info,Opt] = opti_fsolve(fun,x0,opts)
%OPTI_FSOLVE Solve a SNLE using an OPTI SNLE Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_fsolve(fun,x0) solves the nonlinear 
%   equation problem F(x) = 0 where fun is the nonlinear equation vector
%   to be solved, starting at x0.
%
%   [x,fval,exitflag,info] = opti_fsolve(fun,x0,opts) allows the user to 
%   specify optiset options. This includes specifying a solver via the 
%   'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_fsolve(fun,...) returns the internally built 
%   OPTI object.

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)


% Handle missing arguments
if nargin < 3, opts = optiset; end 
if nargin < 2, error('You must supply at least 2 arguments to opti_fsolve'); end

%Sort out Fun + Grad
[f,g] = detGrad(fun,x0);

%Build OPTI Object
if(length(x0) == 1) %scalar problem
    Opt = opti('fun',f,'grad',g,'x0',x0,'probtype','snle','options',opts); 
else
    Opt = opti('fun',f,'grad',g,'x0',x0,'options',opts);
end
    
%Solve
[x,fval,exitflag,info] = solve(Opt);
