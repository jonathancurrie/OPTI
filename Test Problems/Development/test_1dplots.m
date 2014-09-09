%% 1D Plot Testing

%% LP
clc
f = -1;
A = 1;
b = 0.5;
lb = 0;
ub = 5;

Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub)
[x,f,e,i] = solve(Opt)
plot(Opt)

%% LP
clc
f = -1;
A = 1;
b = 0.5;
Aeq = 1; beq = 0.3;
lb = 0;
ub = 5;

Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub)
[x,f,e,i] = solve(Opt)
plot(Opt,[-1 1])

%% MILP
clc
f = -1;
A = 1;
b = 0.5;
lb = 0;
ub = 5;
ivar = 1;

Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'ivars',ivar)
[x,f,e,i] = solve(Opt)
plot(Opt)

%% QP
clc
H = 1
f = -1;
A = 1;
b = 0.5;
lb = 0;
ub = 5;

Opt = opti('H',H,'f',f,'ineq',A,b,'bounds',lb,ub)
[x,f,e,i] = solve(Opt)
plot(Opt)

%% QCQP
clc
H = 1
f = -1;
A = 1;
b = 0.5;
lb = 0;
ub = 5;
Q = 1.5; l = 0.2; r = 0.8;

Opt = opti('H',H,'f',f,'bounds',lb,ub,'qc',Q,l,r)
[x,f,e,i] = solve(Opt)
plot(Opt)

%% QCQP2
clc
H = 1
f = -1;
A = 1;
b = 0.5;
lb = 0;
ub = 5;
Q = 1.5; l = 0.2; r = 0.8;

Opt = opti('H',H,'f',f,'bounds',lb,ub,'qcrow',Q,l,r,r)
[x,f,e,i] = solve(Opt)
plot(Opt)

%% SDP
clc
%Objective 
f = 1;
%Linear Constraints
lb = 0;
ub = 10;
%SDP Constraints x1-5 >= 0
C = 5;
A0 = 1;
sdcone = [C A0];
%Setup Options
opts = optiset('solver','csdp','display','iter');
%Solve
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot Solution
plot(Opt)

%% NLS
clc
%Function
fun = @(x,xdata) x(1)*exp(-0.1013*xdata);
grad = @(x,xdata) [ exp(-0.1013.*xdata), x(1).*xdata.*exp(-0.1013.*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'ndec',1,'options',opts)
x0 = [100]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)


%% Bounded NLP
% Objective
 fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;
% Bounds
 lb = -1.5; 
 ub = 1.5;
% Initial Solution
 x0 = 0;
% Options
 opts = optiset('solver','scip');
% Create OPTI Object
 Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)
% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,[lb,ub])

%% Bounded NLP Multisolve
clc
% Objective
 fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;
% Bounds
 lb = -1.5; 
 ub = 1.5;
% Options
 opts = optiset('solver','ipopt','display','off');
% Create OPTI Object
 Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)
% Solve using SCIP
[x,fval,exitflag,info] = multisolve(Opt)
multiplot(Opt)


%% Bounded Non-differentialable Multi-solve
% Objective
 fun = @RiemannND;
% Bounds
 lb = 0.5;
 ub = 2;
% Starting Guess
 x0 = 1;
% Options
 opts = optiset('solver','nomad','display','iter');
% Create OPTI Object
 Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)
% Attempt to Solve
[x,fval,exitflag,info] = multisolve(Opt,x0)

multiplot(Opt,[],1000)

%% NLP
% Objective
 fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;
% Bounds
 lb = -1.5; 
 ub = 1.5;
%NL con
nlcon = @(x) x^2 - 0.5*x;
nlrhs = 0.5;
nle = -1;
% Initial Solution
 x0 = 0;
% Options
 opts = optiset('solver','scip');
% Create OPTI Object
 Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,[lb,ub])

%% NLP2
% Objective
 fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;
% Bounds
 lb = -1.5; 
 ub = 1.5;
%NL con
nlcon = @(x) [x^2 - 0.5*x; 2*x];
nlrhs = [0.5; 1.5];
nle = [-1;0];
% Initial Solution
 x0 = 0;
% Options
 opts = optiset('solver','scip');
% Create OPTI Object
 Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,[lb,ub])


%% MINLP
% Objective
 fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;
% Bounds
 lb = -1.5; 
 ub = 1.5;
%NL con
nlcon = @(x) x^2 - 0.5*x;
nlrhs = 0.5;
nle = -1;
% Initial Solution
 x0 = 0;
% Options
 opts = optiset('solver','scip');
% Create OPTI Object
 Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ivars',1,'bounds',lb,ub,'options',opts)
% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,[lb,ub])


