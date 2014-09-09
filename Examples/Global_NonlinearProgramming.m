%% OPTI Toolbox Global Nonlinear Programming Examples
%
% This file contains a number of Global NLP problems and demonstrates how 
% to solve them using the OPTI Toolbox. You should read and complete
% BasicUsage.m & NonlinearProgramming.m BEFORE running the below examples.
%
%   Copyright (C) 2014 Jonathan Currie (I2C2)

% There is also a page on the Wiki which supplements this example:
web('http://i2c2.aut.ac.nz/Wiki/OPTI/index.php/Probs/GNLP');

%% Determing which Solver to Use
% OPTI Toolbox comes with a number of NLP solvers, thus to determine which
% ones are available on your system you can type:
clc
checkSolver('NLP')

% Note the columns DR and GL. A cross in DR indicates the solver requires
% 1st (and perhaps 2nd) derivatives, while a cross in GL indicates the
% solver can solve Global Optimization problems. For noisy problems unless
% you have exact derivatives, avoid solvers with a cross in DR.

%% Typical Global Optimization Problems
% Global optimization problems result from any of the following circumstances:
%
%       - Objectives containing noise / a stochastic element
%       - Non-convex functions (functions which are not a bowl / hill in 2D)
%       - Objectives that include periodic functions (sin, cos)
%       - Parameter estimation of ODEs solved with adaptive step integrators
%       - Any problem that contains multiple local minima.
%
% An extreme example from Wolfram is shown below:
clc
%Objective
fun = @(x) norm([x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
          x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))]);
lb = [-4;-4]; ub = [4;4];
x0 = [-4;-4];

%Plot      
n = 1e2;
x = linspace(-4,4,n); y = linspace(-4,4,n); Z = zeros(n,n);
for i = 1:n
    for j = 1:n
        Z(j,i) = fun([x(i),y(j)]);
    end
end
surfc(x,y,Z)
colormap summer; shading interp; lighting phong; view(-38,58);
xlabel('x1'); ylabel('x2'); zlabel('obj'); title('Wolfram Global Optimization Problem'); 

%% Example 1 - Basic Setup
% The main difference when solving a Global NLP is that OPTI treats Global
% and Local NLP problems identically, and therefore you will have to specify
% a Global solver. For this example we will build 4 OPTI objects with 4
% different solvers:
clc
%Build OPTI Problem
prob = optiprob('fun',fun,'bounds',lb,ub);
%Choose Global Solver
opts1 = optiset('solver','nomad');
opts2 = optiset('solver','pswarm');
opts3 = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
opts4 = optiset('solver','ipopt','warnings','off');

%Pass to OPTI Constructor for Error Checking + Setup
Opt1 = opti(prob,opts1); 
Opt2 = opti(prob,opts2); 
Opt3 = opti(prob,opts3); 
Opt4 = opti(prob,opts4); 

% Call solve to solve the problem. Check the plot for a comparison of the
% solution points. Note NOMAD and NLOPT are deterministic with the current
% settings, PSWARM includes random elements, and IPOPT is for comparison of
% a local solution.

[x1,fval1] = solve(Opt1,x0);
[x2,fval2] = solve(Opt2,x0);
[x3,fval3] = solve(Opt3,x0);
[x4,fval4] = solve(Opt4,x0);

view(0,90); hold on;
plot3(x0(1),x0(2),10,'kx','markersize',10);
plot3(x1(1),x1(2),10,'ro'); text(x1(1)+0.1,x1(2)+0.1,10,sprintf('NOMAD: %f',fval1));
plot3(x2(1),x2(2),10,'ro'); text(x2(1)+0.1,x2(2)-0.1,10,sprintf('PSWARM: %f',fval2));
plot3(x3(1),x3(2),10,'ro'); text(x3(1)+0.1,x3(2)+0.2,10,sprintf('NLOPT: %f',fval3));
plot3(x4(1),x4(2),10,'ro'); text(x4(1)+0.1,x4(2)+0.1,10,sprintf('IPOPT: %f',fval4));
hold off;

%% Problem 2 - Quartic
% Includes Linear and Nonlinear Constraints
clc
%Problem
fun = @(x) x(1)^4 - 14*x(1)^2 + 24*x(1) - x(2)^2;

%Linear Constraints
A = [-1 1]; b = 8;
%Nonlinear Constraints
nlcon = @(x) (-x(1)^2) - 2*x(1) + x(2);
nlrhs = -2;
nle = -1;
%Bounds + Starting Guess
lb = [-8;0];
ub = [10;10];
x0 = [0;0];

% Solving a constrained global optimization problem. Note linear
% constraints will be converted to nonlinear ones for this solver
opts = optiset('solver','nomad','solverOpts',nomadset('direction_type','lt 2n')); 
Opt = opti('fun',fun,'ineq',A,b,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)

% This may take a few seconds...
[x,fval,exitflag,info] = solve(Opt,x0)  

%% Example 3 - Solving with PSwarm
% Note x0 is on the local minima side of the saddle
clc
%Problem
fun = @(x) -2*x(1)*x(2);

%Constraints
lb = [-0.5;-0.5];
ub = [1;1];
x0 = [-0.3;-0.3]; 

%Plot      
n = 1e2;
x = linspace(-0.5,1,n); y = linspace(-0.5,1,n); Z = zeros(n,n);
for i = 1:n
    for j = 1:n
        Z(j,i) = fun([x(i),y(j)]);
    end
end
surfc(x,y,Z); hold on; plot3(x0(1),x0(2),fun(x0),'r.','markersize',20); hold off;
colormap winter; shading flat; lighting gouraud; view(18,28);
xlabel('x1'); ylabel('x2'); zlabel('obj'); title('Saddle Point Optimization Problem'); 

% PSwarm solves bounded and linearly constrained global problems
Opt = opti('fun',fun,'bounds',lb,ub,'options',optiset('solver','pswarm'))


[x,fval,exitflag,info] = solve(Opt,x0)  
plot(Opt,3)

%% Example 4 - White Box Quartic
% The following problem is the same as problem 2, however this time we are
% going to solve it using a white box solver (SCIP). The SCIP interface
% will parse the following functions into an algebraic description,
% allowing SCIP to find a global solution to this problem.
clc
%Problem
fun = @(x) x(1)^4 - 14*x(1)^2 + 24*x(1) - x(2)^2;

%Linear Constraints
A = [-1 1]; b = 8;
%Nonlinear Constraints
nlcon = @(x) (-x(1)^2) - 2*x(1) + x(2);
nlrhs = -2;
nle = -1;
%Bounds + Starting Guess
lb = [-8;0];
ub = [10;10];
x0 = [0;0];

% SCIP solves a subset of nonlinear and mixed integer problems, provided
% the problem is deterministics and constains a subset of allowable functions.
Opt = opti('fun',fun,'ineq',A,b,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,...
           'options',optiset('solver','scip'))

[x,fval,exitflag,info] = solve(Opt,x0)  
plot(Opt)

%% Example 5 - OPTI's Multi-Start Solver
% For simple, low dimensional problems, OPTI's mulit-start is a naive
% algorithm for searching the problem space for a global solution. See the
% following page for more examples:

web('http://i2c2.aut.ac.nz/Wiki/OPTI/index.php/Advanced/MultiSolve');

clc
Opt = opti(Opt,'solver','ipopt');
[x,fval,exitflag,info] = multisolve(Opt,x0)  
multiplot(Opt)

%% Summary
% While Global Optimization solvers may take longer, many real engineering
% problems result in noisy or non-convex objectives and local solvers will
% often struggle to return a result, or fall into the closest local
% optima. OPTI provides a range of competitive blackbox global optimization
% solvers which can return much better results, as well as a white box
% solver for academic users which guarantees a global solution.
