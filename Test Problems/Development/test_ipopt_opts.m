%% ALL NONLINEAR VERSION
clc
clear
% Objective & Gradient
 fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
 grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];
 
 % Constraints
A = [-1 1; 1 1];
rl = [-Inf;5];
ru = [-1;5];
lb = [0;0]; ub = [4;4];

% Nonlinear Constraint, Jacobian & Structure
 nlcon = @(x) A*x;
 nljac = @(x) sparse(A);
 jacstr = @() sparse(double(A~=0));        

% Starting Guess
 x0 = [2;2];

% Build Function Structure
 funcs.objective = fun;
 funcs.gradient = grad;
 funcs.constraints = nlcon;
 funcs.jacobian = nljac;
 funcs.jacobianstructure = jacstr;
 funcs.hessian = @(x,sigma,lambda)sparse([[sigma*(1200*x(1)^2-400*x(2)+2),0];[-400*sigma*x(1),200*sigma]]);
 funcs.hessianstructure = @() sparse(tril(ones(2)));

% Build Options Structure
 opts.lb = lb;
 opts.ub = ub;
 opts.cl = rl;
 opts.cu = ru;
%  opts.ipopt.linear_solver = 'pardiso';
%  opts.ipopt.hessian_approximation = 'limited-memory';
opts.ipopt.mu_init = 0.5;

% Call IPOPT
[x,output] = ipopt(x0,funcs,opts)



%% Repeat but with OPTI interface
clc
iopts = ipoptset('mu_init',0.5);
opts = optiset('solverOpts',iopts, 'maxiter',15,'display','iter');
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,rl,ru,'nljac',nljac,'nljacstr',jacstr,'H',funcs.hessian,'hstr',funcs.hessianstructure,'bounds',lb,ub,'opts',opts)

[x,f,e,i] = solve(Opt,x0)













