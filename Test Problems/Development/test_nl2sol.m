%% Test Problems for NL2SOL Interface
% J.Currie 2012
clc
clear all

%% NLS2 [0.3578]
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
% %Setup Options
opts = optiset('solver','nl2sol');
% %Build & Solve
x0=[1.0; 1.0; 1.0];
Opt = opti('fun',fun,'ydata',ydata,'ndec',3,'options',opts)
% [x,fval,exitflag,info] = solve(Opt,x0)


%%
clc
grad = @(x) mklJac(fun,x);

opts = [];
opts.display = 2;
opts.maxiter = 50;
lb = [0;0.5;2];
ub = [10;10;10];

[x,f,e,i,ii] = nl2sol(fun,[],x0,ydata,lb,ub,opts)


%% NLS2 [0.3578]
clc
%Function
i = (1:40)';
fun = @(x) exp(-x(1)*i);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
% %Setup Options
opts = optiset('solver','nl2sol');
% %Build & Solve
x0=[1.0];
Opt = opti('fun',fun,'ydata',ydata,'ndec',1,'options',opts)
% [x,fval,exitflag,info] = solve(Opt,x0)


%%
grad = @(x) mklJac(fun,x);

opts = [];
opts.display = 2;
opts.maxiter = 50;
lb = [0];
ub = [10];

[x,f,e,i,ii] = nl2sol(fun,grad,x0,ydata,[],[],opts)

%%
clc
prob = nls_prob(20);
opts = optiset('solver','nl2sol','display','iter');
Opt = opti(prob,opts)

[x,f,e,i] = solve(Opt)


