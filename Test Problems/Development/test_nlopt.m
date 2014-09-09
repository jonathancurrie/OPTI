%% UNO Rosenbrock
clc
%Objective
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [0 1]';

opts = [];
opts.algorithm = nloptSolver('LN_PRAXIS');
opts.maxfeval = 1000;
opts.maxtime = 1;
opts.display = 2;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% UNO Rosenbrock [GRAD]
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.gradient = @(x) mklJac(nlprob.objective,x);
x0 = [0 1]';

opts = [];
opts.algorithm = nloptSolver('LD_LBFGS');
opts.maxfeval = 100;
opts.maxtime = 1;
opts.display = 2;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% BOUNDED Rosenbrock
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.lb = [-5;-5];
nlprob.ub = [1;0.8];
x0 = [0 0]';

opts = [];
opts.algorithm = nloptSolver('LN_PRAXIS');
opts.maxfeval = 1000;
opts.maxtime = 10;
opts.display = 2;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% BOUNDED GRAD Rosenbrock
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.gradient = @(x) mklJac(nlprob.objective,x);
nlprob.lb = [-5;-5];
nlprob.ub = [1;0.8];
x0 = [0 0]';

opts = [];
opts.algorithm = nloptSolver('LD_LBFGS');
opts.maxfeval = 1000;
opts.maxtime = 10;
opts.display = 2;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% NLCON HS71
clc
%Objective
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);         
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [prod(x); sum(x.^2)];         

nlprob =[];
nlprob.objective = obj;
nlprob.nlcon = nlcon;
nlprob.nlrhs = [25;40];
nlprob.nle = [1;0];
nlprob.lb = lb;
nlprob.ub = ub;

opts = [];
opts.algorithm = nloptSolver('LN_COBYLA');
opts.maxfeval = 1000;
opts.maxtime = 2;
opts.display = 2;
nlprob.options = opts;

x0 = zeros(4,1);

[x,f,e,i] = nlopt(nlprob,x0,opts)

%% NLCON GRAD HS71
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [prod(x); sum(x.^2)];
nljac = @(x) [prod(x)./x'; 2*x'];            

nlprob =[];
nlprob.objective = obj;
nlprob.gradient = grad;
nlprob.nlcon = nlcon;
nlprob.nljac = nljac;
nlprob.nlrhs = [25;40];
nlprob.nle = [1;0];
nlprob.lb = lb;
nlprob.ub = ub;

opts = [];
opts.algorithm = nloptSolver('LD_SLSQP');
opts.maxfeval = 1000;
opts.maxtime = 2;
opts.display = 2;
nlprob.options = opts;

x0 = [1 5 5 1]';

[x,f,e,i] = nlopt(nlprob,x0,opts)

%% NLP1 Hock & Schittkowski #71
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','nlopt','warnings','on','display','iter','solverOpts',nloptset('algorithm','LD_SLSQP'));
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)

%% LP [-31.4]
clc
% clear all
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    

nlprob = [];
nlprob.objective = @(x) f'*x;
nlprob.gradient = @(x) f;

nlprob.nlcon = @(x) A*x;
nlprob.nljac = @(x) A;
nlprob.nlrhs = b;
nlprob.nle = [-1;-1;-1];

nlprob.lb = [0;0];
nlprob.ub = [10;10];

opts = [];
opts.algorithm = nloptSolver('LD_SLSQP');
opts.maxfeval = 100;
opts.maxtime = 2;
opts.display = 2;
nlprob.options = opts;

x0 = [0;0];

[x,f,e,i] = nlopt(nlprob,x0)

%% LP1
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
%Build Object
Opt = opti('obj',@(x) f'*x,'grad',@(x) f,'ineq',A,b,'bounds',[0;0],[10;10],'options',optiset('solver','nlopt','solverOpts',nloptset('algorithm','LD_SLSQP')))
%Build & Solve
[x,fval,exitflag,info] = solve(Opt,[0;0])  
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% Auglag NLP
clc
fun = @(x) log(1+x(1)^2) - x(2);
grad = @(x)[[(2*x(1))/(x(1)^2+1)];[-1]];
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2 - 4;
nljac = @(x)[[4*x(1)*(x(1)^2+1),2*x(2)]];
nlrhs = 0;
nle = 0;
x0 = [2;2];

nlprob = [];
nlprob.objective = fun;
nlprob.gradient = grad;

nlprob.nlcon = nlcon;
nlprob.nljac = nljac;
nlprob.nlrhs = nlrhs;
nlprob.nle = nle;

opts = [];
opts.algorithm = nloptSolver('LD_AUGLAG');
opts.maxfeval = 120;
opts.maxtime = 2;
opts.display = 2;

opts.local_optimizer.algorithm = nloptSolver('LD_LBFGS');
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% Error checking
clc
clear
nlprob.objective = @(x) -1;
nlprob.lb = [-1;-1];
nlprob.ub = [1;1];

nlprob.nlcon = @(x) [-1;-1];
nlprob.nlrhs = [-1;-1];
nlprob.nle = [1;1];

nlprob.options.algorithm = 25;
x0 = [1;1];

[x,f,e,i] = nlopt(nlprob,x0)

%% Stop Val Checking
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.lb = [-5;-5];
nlprob.ub = [1;1];
x0 = [0 0]';

opts = [];
opts.algorithm = nloptSolver('LN_PRAXIS');
opts.maxfeval = 1000;
opts.maxtime = 10;
opts.display = 2;
opts.stopval = 1e-6;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% Vector Storage Checking
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.gradient = @(x) mklJac(nlprob.objective,x);
nlprob.lb = [-5;-5];
nlprob.ub = [1;0.8];
x0 = [0 0]';

opts = [];
opts.algorithm = nloptSolver('LD_LBFGS');
opts.maxfeval = 1000;
opts.maxtime = 10;
opts.display = 2;
opts.vector_storage = 1;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% Initial Pop Testing
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.lb = [0;0];
nlprob.ub = [1;1];
x0 = [0 0]';

opts = [];
opts.algorithm = nloptSolver('GN_CRS2_LM');
opts.maxfeval = 1000;
opts.maxtime = 10;
opts.display = 1;
opts.initial_pop = 30;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)

%% NLP1 Hock & Schittkowski #71 with linear mod
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
A = [1 0 0 1]; rl = 2; ru = 2;
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','nlopt','warnings','on','display','iter','solverOpts',nloptset('algorithm','LD_SLSQP'));
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'lin',A,rl,ru,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)


%% BOUNDED Rosenbrock [ESCH Check]
clc
%Objective
nlprob = [];
nlprob.objective = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
nlprob.lb = [-5;-5];
nlprob.ub = [1;0.8];
x0 = [0.2 0.2]';

opts = [];
opts.algorithm = nloptSolver('GN_ESCH');
opts.maxfeval = 10000;
opts.maxtime = 10;
opts.display = 0;
nlprob.options = opts;

[x,f,e,i] = nlopt(nlprob,x0)