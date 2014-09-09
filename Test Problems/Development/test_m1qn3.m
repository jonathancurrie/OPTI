%% UNO Rosenbrock
clc
clear
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','m1qn3','display','iter');
%Build & Solve
Opt = opti('obj',obj,'ndec',2,'options',opts)
x0 = [0 0]';
[x,fval,exitflag,info]= solve(Opt,x0)

%%
% grad = Opt.nlprob.funcs.gradient;
grad = @(x) autoJac(obj,x);
opts.display = 2;
opts.nupdates = 5;
opts.maxtime = 0.5;
opts.iterfun = [];

[x,f,e,i,ii] = m1qn3(obj,grad,x0,opts);

%% 4x Rosenbrock [1,1,1,1]
clc
%Objective
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
grad = @(x) mklJac(fun,x);
x0 = [randn(1,100) 1];
opts.display = 2;
[x,fval,ef,iter] = m1qn3(fun,grad,x0,opts);