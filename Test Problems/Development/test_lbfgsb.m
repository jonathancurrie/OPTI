%% Test L-BFGS-B

%% 4x Rosenbrock [1,1,1,1]
clc
n = 100;
%Objective
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
grad = @(x) mklJac(fun,x);
x0 = zeros(1,n);
lb = zeros(1,n); ub = 5*ones(1,n);

opts = [];
opts.display = 2;
opts.pgtol = 1e-5;
opts.maxtime = 5;
opts.maxiter = 500;
opts.nupdate = 10;
opts.maxfeval = 90;

[~,fval,ef,iter,feval] = lbfgsb(fun,grad,lb,ub,x0,opts)

%% OPTI VERSION
clc
n = 100;
%Objective (Vectorized Rosenbrock)
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
%Bounds
lb = zeros(1,n); ub = 5*ones(1,n);
x0 = zeros(1,n);
%Setup LBFGSB Options
lopts = lbfgsbset('nupdate',20);
%General Options
opts = optiset('solver','lbfgsb','display','iter','maxtime',5,'maxiter',1500,'maxfeval',1000,'solverOpts',lopts)
%Build and Solve
Opt = opti('fun',fun,'x0',x0,'bounds',lb,ub,'options',opts);
[~,fval,ef,info] = solve(Opt)

%% 4x Rosenbrock [1,1,1,1]
clc
n = 100;
%Objective
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
grad = @(x) mklJac(fun,x);
x0 = zeros(1,n);
lb = zeros(1,n); ub = 5*ones(1,n);

opts = [];
opts.display = 2;
opts.pgtol = 1e-5;
opts.maxtime = 5;
opts.maxiter = 500;
opts.nupdate = 10;
opts.iterfun = @optiplotlogfval;

[~,fval,ef,iter,feval] = lbfgsb(fun,grad,lb,ub,x0,opts)

%% OPTI VERSION
clc
lopts = lbfgsbset('nupdate',20);
opts = optiset('solver','lbfgsb','display','iter','maxtime',5,'maxiter',1500,'maxfeval',1000,'iterfun',@optiplotlogfval,'solverOpts',lopts)
Opt = opti('fun',fun,'x0',x0,'bounds',lb,ub,'options',opts);

[~,fval,ef,info] = solve(Opt)