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
