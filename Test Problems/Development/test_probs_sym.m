%% OPTISYM Test Problems
% A collection of test problems I've used for building the toolbox.
clc
%#ok<*ASGLU,*NASGU,*NOPTS>

%% LP1
clc
%Objective & Constraints
fun = @(x) -[6 5]*x;
con = @(x) [1,4; 6,4; 2, -5]*x;
cl = -Inf(3,1);
cu = [16;28;6];    
lb = [0;0]; ub = [10;10];
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% MILP2
clc
%Objective & Constraints
fun = @(x) -[1 2 3 1]*x; 
con = @(x) [-1 1 1 10; 1 -3 1 0;0 1 0 -3.5]*x; 
cl = [-Inf;-Inf;0];
cu = [20;30;0];
lb = [0 0 0 2]';
ub = [40 inf inf 3]';
xtype = 'CCCI';
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu,xtype)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% BILP1
clc
%Objective & Constraints
fun = @(x) -[6 5]*x;
con = @(x) [-3,5; 6,4; 3, -5; -6, -4]*x; 
cl = -Inf(4,1);
cu = [6;9;1;3];  
xtype = 'BB';
%Build Object
Opt = optisym(fun,[],[],[],con,cl,cu,xtype)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info


%% QP3
clc
%Objective & Constraints
fun = @(x) 0.5*x'*[1 -1; -1 2]*x + -[2 6]*x;
con = @(x) [1 1; -1 2; 2 1; 1 1.5]*x;
cl = [-Inf(3,1);2];
cu = [2; 2; 3; 2];
lb=[0;0]; ub=[10;10];
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% QCQ(L)P2 [-5.2]
clc
%Objective & Constraints
fun = @(x) [-2 -2]*x;
con = @(x) [[-1 1; 1 3]*x; x'*[1 0;0 1]*x + [0 -2]*x];
cl = -Inf(3,1);
cu = [2;5;1];
lb = [0;0];
ub = [40;inf];
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% QCLP3 [-1.7394]
clc
fun = @(x) [-2 -2]*x;
con = @(x) [[-1 1; 1 3]*x; x'*[1 0;0 1]*x + [0 2]*x; x'*[1 0;0 1]*x + [2 -2]*x];
cl = -Inf(4,1);
cu = [2;5;1;1];
lb = [0;0];
ub = [40;inf];
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% Indefinite QCQP [-3.55]
clc
%Objective & Constraints
fun = @(x) 0.5*x.'*eye(2)*x + [-2 -2]*x;
con = @(x) [[-1 1; 1 3]*x; x'*[-1 0;0 1]*x + [0 -2]*x];
cl = -Inf(3,1);
cu = [2;5;-0.5];
lb = [0;0];
ub = [40;inf];
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% MIQCQP1 [-2.5429] [non-convex]
clc
%Objective & Constraints
fun = @(x) 0.5*x.'*eye(2)*x + [-2 -2]*x;
con = @(x) [[-1 1; 1 3]*x; x'*[1 0;0 1]*x + [0 -2]*x];
cl = [-Inf(2,1);3.5];
cu = [2;5;5];
lb = [0;0];
ub = [40;inf];
xtype = 'IC';
%Build Object
Opt = optisym(fun,[],lb,ub,con,cl,cu,xtype)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% NLS1 [9.5049]
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]';
%Function
fun = @(x) sum((x(1)*exp(x(2)*xdata)-ydata).^2);
%Build Object
Opt = optisym(fun,[100;-1])
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% NLP1 Hock & Schittkowski #71 [17.014]
clc
%Objective & Gradient
fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
%Constraints
con = @(x) [ prod(x);
               sum(x.^2)];
cl = [25;40];
cu = [Inf;40];
lb = ones(4,1);
ub = 5*ones(4,1);
x0 = [1 5 5 1]';
%Build Object
Opt = optisym(fun,x0,lb,ub,con,cl,cu)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

%% MINLP1 [-2.1251]
clc
%Objective
fun = @(x) -x(1) - x(2) - x(3)^3;
%Constraints
con = @(x) [ (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.);
                x(1) - x(2);
                x(1) + x(3) + x(4)];
cl = -Inf(3,1);
cu = [1/4;0;2];
ub = [1;Inf;Inf;5];
lb = [0;0;0;0];
xtype = 'BCCI';
x0 = [0;0;0;0];
%Build Object
[Opt,SB] = optisym(fun,x0,lb,ub,con,cl,cu,xtype)
%Build & Solve
[x,fval,exitflag,info] = solve(Opt);
fval
info

