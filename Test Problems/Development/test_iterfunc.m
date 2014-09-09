%% Iteration Callback Function Checking

prob = nlp_prob(50);
opts = optiset('solver','matlab','solverOpts',optimset('PlotFcns',@optimplotfval));
Opt = opti(prob,opts)
[x,f,e,i] = solve(Opt)

%% SNLE
clc
%Equations
f = @(x) [2*x(1) - x(2) - exp(-x(1));
          -x(1) + 2*x(2) - exp(-x(2))];
%Setup Options
opts = optiset('solver','hybrj','display','iter','iterfun',@optiplotfval)
%Build & Solve
x0 = [-5;5];
Opt = opti('fun',f,'x0',x0,'options',opts)
[x,fval,ef,stat] = solve(Opt)

%% NLS1
clc
prob = nls_prob(10);
%Setup Options
opts = optiset('solver','matlab','display','iter','iterfun',@optiplotfval)
%Build & Solve
Opt = opti(prob,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,prob.x0)


%% NLP1
clc
prob = nlp_prob(30);
% prob.fun = @(x) prob.f*x; prob.x0 = zeros(size(prob.f')); prob.f = []; 
%Setup Options
opts = optiset('solver','ipopt','display','iter','iterfun',@optiplotfval)
%Build & Solve
Opt = opti(prob,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,prob.x0)

%% NLP1
clc
prob = nlp_prob(30);
% prob.fun = @(x) prob.f*x; prob.x0 = zeros(size(prob.f')); prob.f = []; 
%Setup Options
opts = optiset('solver','ipopt','display','iter','iterfun',@test_ipoptiter)
%Build & Solve
Opt = opti(prob,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,prob.x0)

%% NLP1 [with fixed variables]
clc
prob = nlp_prob(30);
prob.lb(end) = 2; prob.ub(end) = 2;
% prob.fun = @(x) prob.f*x; prob.x0 = zeros(size(prob.f')); prob.f = []; 
%Setup Options
opts = optiset('solver','ipopt','display','iter','iterfun',@test_ipoptiter)
%Build & Solve
Opt = opti(prob,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,prob.x0)

%% BNLP1
clc
clf
fun = @(x) sum(100*(x(1:end-1).^2-x(2:end)).^2+(x(1:end-1)-1).^2);
n = 100;
lb = -5*ones(n,1); ub = 5*ones(n,1);
x0 = zeros(n,1);

%Setup Options
opts = optiset('solver','lbfgsb','display','iter','iterfun',@optiplotlogfval,'maxfeval',1e6,'maxiter',1e4)
%Build & Solve
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)


%% MINLP3 (not going)
clc
%Objective
obj = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
% Constraints
nlcon = @(x) -x(1)^2 + x(2)-0.5;
nlrhs = 0;
nle = 1;      
%Setup Options
opts = optiset('solver','bonmin','display','iter','iterfun',@optiplotfval,'maxfeval',1e6,'maxiter',1e4)
% Solve
Opt = opti('obj',obj,'ndec',2,'nlmix',nlcon,nlrhs,nle,'int',[1 2],'options',opts)
x0 = [4.9;0.1];
[x,fval,ef,info] = solve(Opt,x0)