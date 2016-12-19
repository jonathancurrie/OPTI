%% Testing PSwarm

%% Rosenbrock [x = 1,1, fval = 0]
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2) - x(1)^2)^2;
x0 = [0 0]';

opts = [];
opts.display = 2;

[x,fval,ef,iter,feval] = pswarm(fun,x0,-[10;10],[10;10],[],[],opts)

%% Constrained Rosenbrock [x = .9488,.9, fval = .0026]
clc
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 0]';
lb = [0;0];
ub = [1;0.9];

opts.display = 2;
opts.maxfeval = 1e6;

[x,fval,ef,iter] = pswarm(fun,x0,lb,ub,[],[],opts)

%% NLP2 Hock & Schittkowski #38
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);                             
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Setup Options
opts = optiset('solver','pswarm','maxiter',1e6,'maxfeval',1e6,'maxtime',0.02);
%Build & Solve
Opt = opti('obj',obj,'bounds',lb,ub,'options',opts)
x0 = [-3  -1  -3  -1]';
[x,fval,exitflag,info]= solve(Opt,x0)