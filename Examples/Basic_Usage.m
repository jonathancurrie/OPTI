%% OPTI Toolbox Basic Functionality Example
%
% This file illustrates the basic functionality of the OPTI Toolbox and
% provides a collection of simple examples to get new users familiar with
% the OPTI methods.
%
% Topics covered include:
%
%   - Checking available solvers
%   - Finding a solver for a particular problem type
%   - Building an optiprob problem structure
%   - Building an optiset options structure
%   - Using the opti constructor
%   - Solving an opti object
%   - Validating the solution
%   - Plotting the solution
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

% There is also a page on the Wiki which supplements this example:
web('https://www.inverseproblem.co.nz/OPTI/index.php/GetStart/Basics');

%% Checking solvers available with your OPTI distribution
% Not all solvers are supplied with the OPTI Toolbox, so to check which are
% available you can use the following function:
clc
optiSolver

%% Matching solvers with a problem type
% Use the 'matrix' argument to print a grid of solver vs problem solved:
clc
optiSolver('matrix')

%% Finding all solvers for a given problem type
% Specify the problem type in order to print a list of all solvers setup to
% solve that problem type. The solver at the top of the list is the default
% solver for that type. The following line prints all solvers which are
% setup to solve a Linear Program (LP) together with constraint and 
% sparsity settings (see help optiSolver):
clc
optiSolver('LP')

%% Creating an Optimization Problem
% To create an optimization problem using OPTI there are two methods
% available:
%
%   1) Pass the problem parameters directly to the opti constructor.
%
%   2) Pass the problem parameters to optiprob to generate a problem
%   structure, then pass the structure to the opti constructor.
%
%   Typically we will use (1), however (2) is provided for backwards
%   compatibility, as well as to break up larger problem descriptions.

%% Problem Fields
% Type opti or optiprob with no input arguments to view all available fields 
% and groupings. Check the help documentation to see alternative naming and 
% grouping strategies.
clc
opti;

%% Options Fields
% Type optiset with no input arguments to view all available fields. Check
% the help documentation for more information.
clc
optiset

web('https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/Opts','-new');

%% Example 1
% This is a simple two decision variable Linear Program (LP) which we will
% use for the next couple of examples.
clc
%Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

%Basic Setup
% Problems can be built using the opti constructor. Options from optiset are 
% always optional. If the user does not provide any options, OPTI will
% choose default ones for the problem being solved.
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub)

%Call solve to solve the problem
[x,fval,exitflag,info] = solve(Opt)   

%% Example 2 - Alternative Setup Strategies
% Naming of arguments, as well as pairing is flexible when using optiprob.
% All example below result in an indentical setup.

Opt = opti('c',f,'ineq',A,b,'bounds',lb,ub);    %c = f

% OR
Opt = opti('grad',f,'ineq',A,b,'bounds',lb,ub); %grad = f

% OR
Opt = opti('f',f,'A',A,'b',b,'bounds',lb,ub);   %individual A,b

% OR
Opt = opti('f',f,'ineq',A,b,'lb',lb,'ub',ub);   %individual lb,ub

%% Example 3 - Choosing the Solver
% You can use optiset to select a solver for your problem. This is supplied
% to the opti constructor using optiset, together with the existing
% problem description.
clc
opts = optiset('solver','clp'); 
Opt = opti(Opt,opts)                            
x = solve(Opt)

%% Example 3b - Choosing the Solver via 'options'
% You can also specify options when creating the initial OPTI object. This
% is the typical way to build models.
clc
opts = optiset('solver','ooqp');
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'options',opts)
x = solve(Opt)

%% Example 4 - Validating the Solution
% The solution is stored inside the OPTI object, so you can validate it
% against the constraints.

[ok,msg] = checkSol(Opt)

%% Example 5 - Plotting the Solution
% Several problem types have a default plot command available IF the
% problem contains two variables.

plot(Opt)

%% Example 6 - Calling a Solver Directly
% Most OPTI solvers use standard MATLAB function prototypes (such as 
% linprog or quadprog), so you can skip using the OPTI class all together:
clc
rl = -Inf(size(b)); %note row format for solver CLP
ru = b;

[x,fval,exitflag,info] = opti_clp([],f,A,rl,ru,lb,ub)

%% Example 7 - MATLAB Optimization Toolbox Overloads
% As described in Overload_demo.m, most Optimization Toolbox functions are
% overloaded ('opti_' is placed in front) for new users.
clc
[x,fval,exitflag,info] = opti_linprog(f,A,b,[],[],lb,ub)


