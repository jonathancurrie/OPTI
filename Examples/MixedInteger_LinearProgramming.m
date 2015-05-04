%% OPTI Toolbox Mixed-Integer Linear Programming Examples
%
% This file contains a number of MILP problems and demonstrates how to 
% solve them using the OPTI Toolbox. You should read and complete
% BasicUsage.m & LinearProgramming.m BEFORE running the below examples.
%
%   Copyright (C) 2014 Jonathan Currie (I2C2)

% There is also a page on the Wiki which supplements this example:
web('http://i2c2.aut.ac.nz/Wiki/OPTI/index.php/Probs/MILP');

%% Determing which Solver to Use
% OPTI Toolbox comes with a number of MILP solvers, thus to determine which
% ones are available on your system you can type:
clc
optiSolver('MILP')

%% Example 1
% This is a simple two decision variable MILP which will use for the next 
% few examples
clc
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];
xtype = 'II';               %Integer Variables (I = integer, C = continuous, B = binary)        

% Building an MILP problem is very similar to an LP, except just add the
% 'xtype' argument for integer variables
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'xtype',xtype)

%Call solve to solve the problem
[x,fval,exitflag,info] = solve(Opt)   

%% Example 2 - Alternative Integer Setup
% You can also supply an vector of integer indices indicating the position of
% continuous and integer variables respectively (note no binary variables 
% can be entered this way)

xtype = [1 2];

Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'xtype',xtype); 
solve(Opt);

%% Example 3 - Plotting the Solution
% Several problem types have a default plot command available. Note for 
% MILP plots it will also plot the integer constraints.

plot(Opt)


%% Example 4 - Solving a slightly bigger MILP
%Build the opti problem:
clc

%Problem
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
lb = [0;0;0;2];
ub = [40;inf;inf;3];
xtype = 'CCCI';

Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'xtype',xtype)

[x,fval,exitflag,info] = solve(Opt) 

%% Example 5 - Calling a solver directly
% All OPTI Linear and Quadratic solvers use standard MATLAB function
% prototypes (such as linprog or quadprog), so you can skip using the OPTI
% class all together. Note be sure to check the integer argument form for
% the particular solver by typing e.g. help opti_glpk
clc
[x,fval,exitflag,info] = opti_glpk(f,A,b,Aeq,beq,lb,ub,xtype)


%% Example 6 - Sparse MILPs
% As with LPs, all solvers are setup to directly solve sparse systems, which 
% is the preferred format for most solvers:
clc
% A larger sparse MILP
load sparseMILP1;

opts = optiset('solver','glpk');    %Solve with GLPK
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'xtype',find(xint),'options',opts)
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% Problem 4
% MILP with Special Ordered Sets (SOS)
clc
%Problem
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
      1 0 1 -3 0];
b = [30;30];
lb = zeros(5,1);
ub = [40;1;inf;inf;1];

%SOS type 1
sos_type = '1';
sos_index = [1 2 3 4 5]';
sos_weight = [1 2 3 4 5]';

% Build the problem, specifying the three SOS fields. Note only some MILP
% solvers are setup to solve problems with SOS:
opts = optiset('solver','cbc');
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'sos',sos_type,sos_index,sos_weight,'options',opts)

[x,fval,exitflag,info] = solve(Opt)
