%% OPTI Toolbox Optimization Toolbox Overloads Examples
%
% This file illustrates the simplest use of OPTI toolbox where a user
% familiar with the Matlab Optimization Toolbox can use OPTI overloads.
%
% In this demo we will cover the following functions:
%
%   - linprog       (Linear Programming)
%   - bintprog      (Binary Integer Programming)
%   - mintprog      (Mixed Integer Linear Programming)
%   - quadprog      (Quadratic Programming)
%   - lsqcurvefit   (Nonlinear Least Squares)
%   - fsolve        (System of Nonlinear Equations)
%   - fminunc       (Unconstrained Nonlinear Optimization)
%   - fmincon       (Constrained Nonlinear Optimization)
%
% As well as setting up OPTI solver options which can be used with the
% overloaded functions.
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

% There is also a page on the Wiki which supplements this example:
web('https://www.inverseproblem.co.nz/OPTI/index.php/GetStart/Overloads');

%% Example 1 - opti_linprog()
% Solved using an OPTI LP solver. Note the function prototype is identical
% to that of the Optimization Toolbox function, but with 'opti_' in front
% of it:
clc
%Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

[x,fval,exitflag,info] = opti_linprog(f,A,b,[],[],lb,ub)


%% Example 2 - opti_bintprog()
% Solved using an OPTI BILP/MILP solver.
clc
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];

[x,fval,exitflag,info] = opti_bintprog(f,A,b)


%% Example 3 - opti_minprog()
% Solved using an OPTI MILP solver. Note this function existed before
% intlinprog.
clc

f = [2, 3, 7, 7];
A = [-1, -1, 2, 5;1, -2, -1, -4];
b = [-2, -3]';
lb = zeros(4,1); 
ub = [30 100 20 1]';  
xint = 'CICI';

[x,fval,exitflag,info] = opti_mintprog(f,A,b,[],[],lb,ub,xint)

%% Example 3b - opti_intlinprog()
% Solved as above. Note that integer variables should be given as integer
% indices.
clc

f = [2, 3, 7, 7];
A = [-1, -1, 2, 5;1, -2, -1, -4];
b = [-2, -3]';
lb = zeros(4,1); 
ub = [30 100 20 1]';  
xint = [2 4];

[x,fval,exitflag,info] = opti_intlinprog(f,xint,A,b,[],[],lb,ub)


%% Example 4 - opti_quadprog()
% Solved using an OPTI QP solver.
clc

H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];

[x,fval,exitflag,info] = opti_quadprog(H,f,A,b)

%% Example 5 - opti_lsqcurvefit()
% Solved using an OPTI NLS solver.
clc

fun = @(x,xdata) x(1)*exp(x(2)*xdata);
x0 = [100; -1];
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

[x,fval,exitflag,info] = opti_lsqcurvefit(fun,x0,xdata,ydata)


%% Example 6 - opti_fsolve()
% Solved using an OPTI SNLE solver.
clc

fun = @(x) [2*x(1) - x(2) - exp(-x(1));
             -x(1) + 2*x(2) - exp(-x(2))];
x0 = [-5;5];

[x,fval,exitflag,info] = opti_fsolve(fun,x0)

%% Example 7 - opti_fminunc()
% Solved using an OPTI UNO solver.
clc

fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';

[x,fval,exitflag,info] = opti_fminunc(fun,x0)


%% Example 8 - opti_fmincon()
% Solved using an OPTI NLP solver.
clc

fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;    
lb = [-1.5;-3];
ub = [4;3];
x0 = [0;0];

[x,fval,exitflag,info] = opti_fmincon(fun,x0,[],[],[],[],lb,ub)

%% Adding Options
% You may also add OPTI options to the overloaded function. For example to
% choose an alternative NLP solver other than the default, examine the
% other solvers available:
clc
optiSolver('NLP')

% Then select it via optiset:

opts = optiset('solver','lbfgsb')

%% Example 9 - NLP with L-BFGS-B
% Now the previous bounded problem can be solved with the L-BFGS-B solver
% by passing 'opts' to the routine:
clc
[x,fval,exitflag,info] = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)

%% Example 10 - Displaying solver information
% You can also display iteration by iteration information from some
% solvers by setting the 'display' field:
clc
opts = optiset(opts,'display','iter');

x = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)

