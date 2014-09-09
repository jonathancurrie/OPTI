clc

% Objective (Nonlinear Equations) Function
 nleq = @(x) [ 2*x(1) - x(2) - exp(-x(1));
             -x(1) + 2*x(2) - exp(-x(2))];

nlJac = @(x) sparse([[exp(-x(1))+2,-1];[-1,exp(-x(2))+2]]);         
         
         
% Starting Guess
 x0 = [-5;5];

% Create OPTI Object
 Opt = opti('nleq',nleq,'nlJac',nlJac,'x0',x0)

% Solve the SNLE problem
[x,fval,exitflag,info] = solve(Opt)