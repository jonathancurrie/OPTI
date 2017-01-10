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

% Build Options Structure
 opts.lb = lb;
 opts.ub = ub;
 opts.cl = rl;
 opts.cu = ru;
 opts.ipopt.hessian_approximation = 'limited-memory';

% Call IPOPT
[x,output] = bonmin(x0,funcs,opts)


%% ALL LINEAR VERSION
clc
clear all
% Objective & Gradient
 fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
 grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];
 
 % Constraints
A = [-1 1; 1 1];
rl = [-Inf;5];
ru = [-1;5];
lb = [0;0]; ub = [4;4];

% Starting Guess
 x0 = [2;2];

% Build Function Structure
 funcs.objective = fun;
 funcs.gradient = grad;

% Build Options Structure
 opts.lb = lb;
 opts.ub = ub;
 opts.A = sparse(A);
 opts.rl = rl;
 opts.ru = ru;
 opts.ipopt.hessian_approximation = 'limited-memory';
 opts.display = 2;

% Call IPOPT
[x,output] = bonmin(x0,funcs,opts)

%% NL INEQ, LIN EQ
clc
clear all
% Objective & Gradient
 fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
 grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];
 
 % Constraints
Aeq = [1 1];
rl = 5;
ru = 5;
lb = [0;0]; ub = [4;4];

% Nonlinear Constraint, Jacobian & Structure
 A = [-1 1];
 nlcon = @(x) [A]*x;
 nljac = @(x) sparse([A]);
 jacstr = @() sparse(double([A]~=0)); 
 cl = -Inf;
 cu = -1;

% Starting Guess
 x0 = [2;2];

% Build Function Structure
 funcs.objective = fun;
 funcs.gradient = grad;
 funcs.constraints = nlcon;
 funcs.jacobian = nljac;
 funcs.jacobianstructure = jacstr;

% Build Options Structure
 opts.lb = lb;
 opts.ub = ub;
 opts.A = sparse(Aeq);
 opts.rl = rl;
 opts.ru = ru;
 opts.cl = cl;
 opts.cu = cu;
 opts.ipopt.hessian_approximation = 'limited-memory';
  opts.display = 2;

% Call IPOPT
% for i = 1:1000
    [x,output] = bonmin(x0,funcs,opts)
% end

%% NL EQ, LIN INEQ
clc
clear all
% Objective & Gradient
 fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
 grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2,200*x(2)-200*x(1)^2];
 
 % Constraints
A = [-1 1];
rl = -Inf;
ru = -1;
lb = [0;0]; ub = [4;4];

% Nonlinear Constraint, Jacobian & Structure
 Aeq = [1 1];
 nlcon = @(x) [Aeq]*x;
 nljac = @(x) sparse([Aeq]);
 jacstr = @() sparse(double([Aeq]~=0)); 
 cl = 5;
 cu = 5;

% Starting Guess
 x0 = [2;2];

% Build Function Structure
 funcs.objective = fun;
 funcs.gradient = grad;
 funcs.constraints = nlcon;
 funcs.jacobian = nljac;
 funcs.jacobianstructure = jacstr;

% Build Options Structure
 opts.lb = lb;
 opts.ub = ub;
 opts.A = sparse(A);
 opts.rl = rl;
 opts.ru = ru;
 opts.cl = cl;
 opts.cu = cu;
 opts.ipopt.hessian_approximation = 'limited-memory';
  opts.display = 1;

% Call IPOPT
% for i = 1:1000
    [x,output] = bonmin(x0,funcs,opts)
% end

%%
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
A = [-1 1]; 
b = -1;
Aeq = [1 1]; 
beq = 5; 
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
opts = optiset('solver','bonmin','display','iter');
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ineq',A,b,'eq',Aeq,beq,'ivars',2,'options',opts)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
% plot(Opt,[],1)

%%
% Objective & Gradient
fun = @(x) -sum(x);
grad = @(x) [-1 -1 -1 -1];

% Linear constraints
n = 4;
A = randn(n,4);
rl = zeros(n,1);
ru = zeros(n,1);

% Nonlinear Constraint, Jacobian & Structure 
nlcon = @(x) x(:).^4; 
nljac = @(x) sparse(diag(4*x(:).^3)); 
jacstr = @() speye(4); 
cl = [-inf;-inf;-inf;-inf]; 
cu = [1;1;1;1];

% Starting Guess
x0 = [0;0;0;0];
% Build Function Structure
funcs.objective = fun;
funcs.gradient = grad;
funcs.constraints = nlcon;
funcs.jacobian = nljac;
funcs.jacobianstructure = jacstr;

% Build Options Structure
lb = [-10;-10;-10;-10]; 
ub = [40;40;40;40]; 
opts.lb = lb; 
opts.ub = ub; 
opts.rl = rl; 
opts.ru = ru; 
opts.cl = cl; 
opts.cu = cu; 
opts.A = sparse(A); 
opts.ipopt.hessian_approximation = 'limited-memory';
opts.display = 2;

% Call IPOPT
[x,output] = bonmin(x0,funcs,opts)



%% SOS Testing - 1
clc

%Problem Data
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;100;100;1];

sos = [];
sos.type = '1';
sos.index = 1:5;
sos.weight = sos.index;

% Objective & Gradient
fun = @(x) f'*x;
grad = @(x) f';

% Linear constraints
rl = -100*ones(size(b));
ru = b;

% Starting Guess
x0 =zeros(size(f));
% Build Function Structure
funcs = [];
funcs.objective = fun;
funcs.gradient = grad;

% Build Options Structure 
opts = [];
opts.lb = zeros(size(ub)); 
opts.ub = ub; 
opts.rl = rl; 
opts.ru = ru; 
opts.A = sparse(A); 
opts.sos = sos;
opts.ipopt.hessian_approximation = 'limited-memory';
opts.display = 2;

% Call BONMIN
[x,output] = bonmin(x0,funcs,opts)

%Check it
Opt = opti('f',f,'lin',A,rl,ru,'bounds',opts.lb,opts.ub,'sos',sos);
[x,fval] = solve(Opt)


%% SOS Testing - 2
clc

%Problem Data
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];

sos = [];
sos.type = '2';
sos.index = [3:5]';
sos.weight = [3:5]';

% Objective & Gradient
fun = @(x) f'*x;
grad = @(x) f';

% Linear constraints
rl = -Inf(size(b));
ru = b;

% Starting Guess
x0 =zeros(size(f));
% Build Function Structure
funcs = [];
funcs.objective = fun;
funcs.gradient = grad;

% Build Options Structure 
opts = [];
opts.lb = -Inf(size(ub)); 
opts.ub = ub; 
opts.rl = rl; 
opts.ru = ru; 
opts.A = sparse(A); 
opts.sos = sos;
opts.ipopt.hessian_approximation = 'limited-memory';
opts.display = 2;

% Call BONMIN
[x,output] = bonmin(x0,funcs,opts)

%Check it
Opt = opti('f',f,'lin',A,rl,ru,'bounds',opts.lb,opts.ub,'sos',sos);
[x,fval] = solve(Opt)


%% SOS Testing - 12
clc

%Problem Data
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];

sos = [];
sos.type = '12';
sos.index = {[1:2]' [3:5]'};
sos.weight = {[1:2]' [3:5]'};

% Objective & Gradient
fun = @(x) f'*x;
grad = @(x) f';

% Linear constraints
rl = -Inf(size(b));
ru = b;

% Starting Guess
x0 =zeros(size(f));
% Build Function Structure
funcs = [];
funcs.objective = fun;
funcs.gradient = grad;

% Build Options Structure 
opts = [];
opts.lb = -Inf(size(ub)); 
opts.ub = ub; 
opts.rl = rl; 
opts.ru = ru; 
opts.A = sparse(A); 
opts.sos = sos;
opts.ipopt.hessian_approximation = 'limited-memory';
opts.display = 2;

% Call BONMIN
[x,output] = bonmin(x0,funcs,opts)

%Check it
Opt = opti('f',f,'lin',A,rl,ru,'bounds',opts.lb,opts.ub,'sos',sos);
[x,fval] = solve(Opt)


%% OPTI call version of SOS 12
clc
opts = optiset('display','iter','solver','bonmin','warnings','all');
Opt = opti('fun',fun,'lin',A,rl,ru,'bounds',zeros(size(ub)),ub,'sos',sos,'opts',opts);
[x,fval] = solve(Opt,x0)
