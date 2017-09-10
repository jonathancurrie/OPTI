%% NLS - MEX Interface
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
% Fitting function
fun = @(x) x(1)*exp(x(2)*xdata);
grad = @(x) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
x0 = [100; -1]; % Starting guess

% Build gsl problem struct
prob = [];
prob.fun        = fun;
prob.grad       = grad;
prob.x0         = x0;
prob.ydata      = ydata;
prob.probType   = 'nls';
prob.options    = [];
prob.options.display = 2;
prob.options.scalingMethod = 'more';
prob.options.linearSolver = 'cholesky';
prob.options.trustRegionSolver = 'dogleg';

[x,f,e,i] = gsl(prob)





%% NLS1
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

%Setup Options
opts = optiset('solver','gsl','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess

[x,f,e,i] = solve(Opt,x0)


%% NLS2
clc
%Function
xdata = (1:40)';
fun = @(x,xdata) x(1)*exp(-x(2)*xdata) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Setup Options
opts = optiset('solver','gsl','display','iter','solverOpts',mkltrnlsset('tolTR',1e-9));
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'data',xdata,ydata,'ndec',3,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt)

%% NLS Example Prob
clc
%Get Problem
prob = nls_prob(2);
opts = optiset('display','iter','solver','gsl','maxfeval',1e5,'maxiter',1e4);
%Build OPTI Object
Opt = opti(prob,opts);
%Solve
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt,[],1)

%% UNO Rosenbrock
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','gsl','display','iter');
%Build & Solve
Opt = opti('obj',obj,'ndec',2,'options',opts)
x0 = [0 0]';
[x,fval,exitflag,info]= solve(Opt,x0)
%Plot
plot(Opt,[],1)

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
opts = optiset('solver','ipopt','warnings','all','display','iter','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)
info.Lambda

%% NLP2 Hock & Schittkowski #38
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)]';
%Hessian (must not be tril for matlab!) & Hessian Structure       
hess = @(x) sparse(  [ 1200*x(1)^2-400*x(2)+2  -400*x(1)       0                          0
                       -400*x(1)               220.2           0                          19.8
                        0                      0               1080*x(3)^2- 360*x(4) + 2  -360*x(3)
                        0                      19.8            -360*x(3)                   200.2 ]);  
Hstr = @() sparse([ 1  1  0  0 
                    1  1  0  1
                    0  0  1  1
                    0  1  1  1 ]);                               
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Setup Options
opts = optiset('solver','ipopt','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub,'options',opts)
x0 = [-3  -1  -3  -1]';
[x,fval,exitflag,info]= solve(Opt,x0)

%% NLP3 Hock & Schittkowski #51
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ]';
%Nonlinear Constraints & Jacobian Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) sparse([ 1  3  0  0  0;
	                  0  0  1  1 -2;
	                  0  1  0  0 -1 ]);
nljacstr = @() sparse([1 1 0 0 0;
                       0 0 1 1 1;
                       0 1 0 0 1]);
nlrhs = [4 0 0]';
nle = [0 0 0]';
%Setup Options
opts = optiset('solver','nlopt','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr,'options',opts)
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,exitflag,info]= solve(Opt,x0)

% %% NLP4 Hock & Schittkowski #104 (not working)
% clc
% %Objective
% obj = @(x) 0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
% %Constraints
% nlcon = @(x) [ 1 - 0.588*x(5)*x(7) - 0.1*x(1);
%                1 - 0.588*x(6)*x(8) - 0.1*x(1) - 0.1*x(2);
%                1 - 4*x(3)*x(5)^(-1) - 2*x(3)^(-0.71)*x(5)^(-1) - 0.588*x(3)^(-1.3)*x(7);
%                1 - 4*x(4)*x(6)^(-1) - 2*x(4)^(-0.71)*x(6)^(-1) - 0.588*x(4)^(-1.3)*x(8);
%                0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
%                0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2)];
% nlrhs = [ 0 0 0 0 1 4.2]';
% nle = [1 1 1 1 1 -1]';
% lb = 0.1*ones(8,1);
% ub = 10*ones(8,1);
% %Setup Options
% opts = optiset('solver','ipopt','display','iter');
% %Build & Solve
% x0 = [6,3,.4,.2,6,6,1,.5];
% Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
% [x,fval,exitflag,info]= solve(Opt,x0)

%% NLP HS Interface
clc
%Get HS problem
no = 15;
[prob,sol,fmin] = nlp_HS(no);
%Build & Solve
Opt = opti(prob)
[x,fval,exitflag,info]= solve(Opt)

acc = norm(fmin-fval)

plot(Opt,[],1) %problems 1-24 can be plotted

%% NLP Riemann Function
clc
%Function in file
fun = @RiemannND;
%Bounds
lb = 3; ub = 5;
x0 = 2.5;

Opt = opti('fun',fun,'x0',x0,'bounds',lb,ub,'options',optiset('solver','nomad'))
[x,fval,exitflag,info]= solve(Opt)

plot(Opt,[2 6],[],1000)