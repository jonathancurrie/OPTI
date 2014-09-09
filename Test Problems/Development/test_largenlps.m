%% Large NLP
clc
clear
%Number of decision vars AND number of linear constraints AND number of
%nonlinear constraints
n = 2000;

% Objective & Gradient
fun = @(x) -sum(x);
grad = @(x) -ones(n,1);

% Linear constraints
A = randn(n,n);
rl = zeros(n,1);
ru = zeros(n,1);
%%
% Nonlinear Constraint, Jacobian & Structure 
nlcon = @(x) x(:).^4; 
nljac = @(x) sparse(diag(4*x(:).^3)); 
jacstr = @() speye(n); 
hessian = @(x,sigma,lambda) sparse(diag(12*x(:).^2));
hessianstructure = @() speye(n);
cl = -Inf(n,1); 
cu = ones(n,1);

% Starting Guess
x0 = 1*ones(n,1);
% Build Function Structure
funcs.objective = fun;
funcs.gradient = grad;
funcs.constraints = nlcon;
funcs.jacobian = nljac;
funcs.jacobianstructure = jacstr;
% funcs.hessian = hessian;
% funcs.hessianstructure = hessianstructure;

% Build Options Structure
lb = -10*ones(n,1); 
ub = 40*ones(n,1); 
opts = [];
opts.lb = lb; 
opts.ub = ub; 
opts.rl = rl; 
opts.ru = ru; 
opts.cl = cl; 
opts.cu = cu; 
opts.A = sparse(A); 
opts.ipopt.hessian_approximation = 'limited-memory';
opts.ipopt.linear_solver = 'pardiso';
% Call IPOPT
[~,output] = ipopt(x0,funcs,opts);


%% Modify to call BONMIN
clc
opts.display = 1;
opts.var_type = zeros(n,1); opts.var_type(1:10) = 1;
[~,output] = bonmin(x0,funcs,opts);




