%% OPTI Toolbox Quadratic Programming Examples
%
% This file contains a number of QP problems and demonstrates how to 
% solve them using the OPTI Toolbox. You should read and complete
% BasicUsage.m & LinearProgramming.m BEFORE running the below examples.
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

% There is also a page on the Wiki which supplements this example:
web('https://www.inverseproblem.co.nz/OPTI/index.php/Probs/QP');

%% Determing which Solver to Use
% OPTI Toolbox comes with a number of QP solvers, thus to determine which
% ones are available on your system you can type:
clc
optiSolver('QP')

%% Example 1
% This is a simple two decision variable QP which will use for the next 
% few examples
clc
%Problem
H = [1 -1; -1 2];           %Objective Function (min 1/2x'Hx + f'x)
f = -[2 6]';
A = [1 1; -1 2; 2 1];       %Linear Inequality Constraints (Ax <= b)
b = [2; 2; 3];
lb = [0;0];                 %Bounds on x (lb <= x <= ub)    


% Building an QP problem is very similar to an LP, except just add the
% H argument for the problem quadratic objective
Opt = opti('H',H,'f',f,'ineq',A,b,'lb',lb)


%Call solve to solve the problem
[x,fval,exitflag,info] = solve(Opt)   

%% Example 2 - Alternative Setup Strategies
% Naming of arguments, as well as pairing is flexible when using optiprob
clc
Opt = opti('hess',H,'f',f,'ineq',A,b,'lb',lb) %hess = H

% OR
Opt = opti('H',H,'grad',f,'ineq',A,b,'lb',lb) %grad = f

%OR
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb)

%% Example 3 - Plotting the Solution
% Several problem types have a default plot command available IF the
% problem contains two variables.

solve(Opt);
plot(Opt)

%% Example 4 - Solving an unconstrained QP
% An unconstrained QP (effectively solving -H\f)
clc
%Problem
H = [1 -1; -1 2];
f = -[2 6]';

Opt = opti('hess',H,'grad',f)

[x,fval,exitflag,info] = solve(Opt)

%% Example 5 - Indefinite QP
% OPTI QP solvers only solve Positive Definite problems. However you can
% have a go with indefinite problems, but you may only get a local
% solution.
clc
%Problem
H = [0 -2; -2 0];           %Objective Function (min 1/2x'Hx + f'x)
f = [0 0]';
lb = [-0.5;-0.5];           %Bounds on x (lb <= x <= ub)   
ub = [1;1];

% Options
opts = optiset('solver','ooqp');

% Create OPTI Object
Opt = opti('qp',H,f,'bounds',lb,ub,'options',opts)

% Solve the Indefinite QP
[x,f] = solve(Opt)

plot(Opt,2)
