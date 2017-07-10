%% OPTI Test Problems
% A collection of test problems I've used for building the toolbox.
clc
%Run the following line to see what solvers are available to OPTI:
optiSolver();

%Alternatively you can find the best available solver for a given problem 
%type:
lp = optiSolver('best_lp')

%Or see all available solvers for a given problem type
optiSolver('lp')

%#ok<*ASGLU,*NASGU,*NOPTS>

%% LP1
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
%Build Object
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10])
%Build & Solve
[x,fval,exitflag,info] = solve(Opt)  
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1 1],40,'bounds',[0 0 0]',[40 inf inf]','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% LP3
clc
%Objective & Constraints
f = [8,1]';
A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
b = [-4,2,21,39,3,0,0]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'options',opts)
x0 = [9,6]';
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt,10)

%% SLE
clc
clear
%Problem
n = 800;
A = sprandn(n,n,n*n/10);
b = 1:n; b = b';
%Setup Options
opts = optiset('solver','mumps');
%Build & Solve
Opt = opti('sle',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt);
info

%% SLE2
clc
clear
%Problem
n = 2;
A = randn(n,n);
b = randn(n,1);
%Setup Options
opts = optiset('solver','mumps');
%Build & Solve
Opt = opti('sle',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];  
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP2
clc
%Objective & Constraints
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',[0 0 0 2]',[40 inf inf 3]','int','CCCI','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP3 Infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = [4,1; 1.5,1; -2, 1; -0.2, 1]; 
b = [5.3;4.5;-2.5; 1.5];  
%Setup Options
opts = optiset('solver','lp_solve');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% MILP4
clc
%Objective & Constraints
f = -[-1, 2]';
A = [2, 1;-4, 4];
b = [5, 5]';
e = -[1, 1];
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('grad',f,'mix',A,b,e,'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% MILP5
clc
%Objective & Constraints
f = [3, -7, -12]';
A = [-3, 6, 8;6, -3, 7;-6, 3, 3];
b = [12, 8, 5]';
e = [-1, -1, -1];
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('grad',f,'mix',A,b,e,'int','III','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% MILP6
clc
%Objective + Constraints
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];
%Special Ordered Sets
sos = '12';
sosind = {(1:2) (3:5)};
soswt = {(1:2) (3:5)};
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('f',f,'ineq',A,b,'ub',ub,'sos',sos,sosind,soswt,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% BILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BB','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt,3)

%% BILP2
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
b = [9; 1; 0; 0];
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BBBB','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% QP1
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','ooqp');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% QP2
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
%Setup Options
opts = optiset('display','iter','solver','auto');
%Build & Solve
Opt = opti('qp',H,f,'ineq',A,b,'lb',[0;0],'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% QP3
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
Aeq = [1 1.5];
beq = 2;
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('qp',H,f,'eq',Aeq,beq,'ineq',A,b,'bounds',[0;0],[10;10],'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% QP4 unconstrained
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
%Build & Solve
Opt = opti('hess',H,'grad',f)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% Nonconvex QP
clc
%Objective & Constraints
H = sparse([0 -1; -1 0]);
f = [0;0];
%Constraints
lb = [-0.5;-0.5];
ub = [1;1];
%Build & Solve
Opt = opti('qp',H,f,'bounds',lb,ub,'options',optiset('solver','auto','display','iter'))
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% QCQP1 [-0.4004]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [20;30];
Q = eye(3);
l = [2;2;2];
r = 0.5;
lb = [0;0;0];
ub = [40;inf;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('display','iter','solver','ipopt'))
[x,fval,exitflag,info] = solve(Opt)

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% QCQ(L)P2 [-5.2]
clc
%Objective & Constraints
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% QCQP3 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1]};
l = [0 2;2 -2];
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% QCQP3a [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
%Alternate QC form
Q = {[1 0;0 1];[1 0;0 1]};
l = {[0;2];[2 -2]};
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% QCQP4 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1];[1 0;0 1]};
l = [0 2 -2;2 -2 -1];
r = [1 1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% Positive Semidefinite QCQP [-3.0834]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;3.8];
Q = [0 0;0 10];
l = [0;-2];
r = 3;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('solver','auto'))
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt,3)

%% Indefinite QCQP [-3.55]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [-1 0;0 1];
l = [0;-2];
r = -0.5;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('solver','scip'))
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% MIQP1 [-11.5]
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
b = [3; 3.5];
ivars = 1;
%Setup Options
opts = optiset('display','iter','solver','bonmin');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int',ivars,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,4)

%% MIQP2 [-2.75]
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'int','CIC','options',opts)
try
[x,fval,exitflag,info] = solve(Opt)
catch
end

%% MIQCQP1 [-2.5429] [non-convex]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
qrl = 3.5;
qru = 5;
lb = [0;0];
ub = [40;inf];
xtype = 'IC';
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,qrl,qru,'bounds',lb,ub,'xtype',xtype)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% MIQCQP2 [-0.0942]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [20;30];
Q = eye(3);
l = [0;0;0];
r = 1;
lb = [0;0;0];
ub = [40;inf;inf];
int = 'CCI';
%Setup Options
opts = optiset('solver','auto','display','iter');
% Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'int',int,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

[ok,msg] = checkSol(Opt)

%% SDP1 [1.4142]
%Note OPTI accepts Standard Primal Form [sum(xi*Ai) - C >= 0]
clc
%Objective
f = 1;
%SDP Constraints [x sqrt(2); sqrt(2) x] >= 0
C = -[0 sqrt(2); sqrt(2) 0];
A = eye(2); 
sdcone = [C A];
%Setup Options
opts = optiset('solver','csdp','display','iter');
% Solve
Opt = opti('f',f,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

[ok,msg] = checkSol(Opt)
plot(Opt)

%% SDP2 [4]
clc
%Objective 
f = [1;1];
%Linear Constraints
lb = [0;0];
ub = [10;10];
%SDP Constraints [x1 2; 2 x2] >= 0
C = -[0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
sdcone = {[C A0 A1]};
%Setup Options
opts = optiset('solver','dsdp','display','iter');
%Solve
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot Solution
plot(Opt)

%% SDP3 [1.2] (Note assume matrices are symmetric triu to avoid x5)
clc
clear
%Objective
f = [1 0 0 0]';
%SDP Constraint1 [x2 x3; x3 x4] <= x1*eye(2)
C = zeros(2);
A0 = eye(2);
A1 = -[1 0; 0 0];
A2 = -[0 1; 1 0];
A3 = -[0 0; 0 1];
sdcone{1} = [C A0 A1 A2 A3];
%SDP Constraint2 [x2 x3; x3 x4] >= [1 0.2; 0.2 1]
C  = [1 0.2; 0.2 1];
A0 = zeros(2);
A1 = [1 0; 0 0];
A2 = [0 1; 1 0];
A3 = [0 0; 0 1];
sdcone{2} = [C A0 A1 A2 A3];
%Setup Options
opts = optiset('solver','dsdp','display','iter');
% Solve
Opt = opti('f',f,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
[ok,msg] = checkSol(Opt)
plot(Opt,1)

%% SDP4 [10.1787]
clc
clear
%Objective
f = [1 1 1];
%Linear Constraints
lb = [10;0;0];
ub = [1000;1000;1000];
%SDP Constraints [x1 1 2; 1 x2 3; 2 3 100] >= 0
C = -[0 1 2; 1 0 3; 2 3 100];
A0 = [1 0 0; 0 0 0; 0 0 0];
A1 = [0 0 0; 0 1 0; 0 0 0];
A2 = zeros(3);
sdcone = [C A0 A1 A2];
%Setup Options
dopts = dsdpset('ptol',1e-8,'rho',3,'zbar',0);
opts = optiset('solver','dsdp','display','iter','solverOpts',dopts);
% Solve
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% SDP5 [Sparse Column Format] [4]
clc
%Objective 
f = [1;1];
%Linear Constraints
lb = [0;0];
ub = [10;10];
%SDP Constraints [x1 2; 2 x2] >= 0
C = -[0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
% sdcone = [C A0 A1];
sdcone = sparse([C(:) A0(:) A1(:)]);
%Setup Options
opts = optiset('solver','csdp','display','iter');
% Solve
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% SDP6 [SeDuMi Structure Format] [sqrt(2)]
clc
%Objective 
b = -1;
%SDP Constraints [x sqrt(2); sqrt(2) x] >= 0
C = -[0 sqrt(2); sqrt(2) 0];
A = -eye(2); 
K = struct('s',2);
%Setup Options
opts = optiset('solver','csdp','display','iter');
% Solve
Opt = opti('sedumi',A(:),b,C(:),K,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% SDP7 [SeDuMi Structure with Linear Con] [4.1667]
clc
%Objective 
b = -[1;1];
%Linear Constraints ([0;0] <= x <= [10;1.5])
c = [0;0;10;1.5];
A = [-1 0; 0 -1; 1 0; 0 1];
%SDP Constraints [x1 2; 2 x2] >= 0
C = -[0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
%Concantenate
A = [A; -A0(:) -A1(:)];
c = [c; -C(:)];
K = struct('l',4,'s',2);
sdcone = [];
sdcone.At = A;
sdcone.b = b;
sdcone.c = c;
sdcone.K = K;
%Setup Options
opts = optiset('solver','csdp','display','iter');
% Solve
Opt = opti('sedumi',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% SNLE1
clc
%Equations
nleq = @(x) x^2 - 4;
%Build & Solve
x0 = 1;
Opt = opti('nleq',nleq,'x0',x0,'options',optiset('solver','auto','display','iter'))
[x,fval,ef,stat] = solve(Opt)
plot(Opt)

%% SNLE2
clc
%Equations
nleq = @(x) [2*x(1) - x(2) - exp(-x(1));
          -x(1) + 2*x(2) - exp(-x(2))];
%Build & Solve
x0 = [-5;5];
Opt = opti('nleq',nleq,'x0',x0,'options',optiset('solver','auto','display','iter'))
[x,fval,ef,stat] = solve(Opt)
plot(Opt,[],1)

%% SNLE3
clc
%Equations
nleq = @(x) [10*(x(2) - x(1)^2)
             sqrt(90)*(x(4) - x(3)^2)
             sqrt(10)*(x(2) + x(4) - 2)
             (1/sqrt(10))*(x(2) - x(4))];
%Build & Solve
x0 = [-30;-10;-30;-10];
Opt = opti('nleq',nleq,'x0',x0,'options',optiset('solver','auto','display','iter'))
[x,fval,ef,stat] = solve(Opt)
plot(Opt,3,1)

%% SCNLE1
clc
%Equations
nleq = @(x) [2*x(1) - x(2) - exp(-x(1));
          -x(1) + 2*x(2) - exp(-x(2))];
%Build & Solve
x0 = [-5;5];
Opt = opti('nleq',nleq,'x0',x0,'bounds',[0.6;0],[1;1],'options',optiset('solver','auto','display','iter'))
[x,fval,ef,stat] = solve(Opt)
plot(Opt,1)

%% NLS1
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

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
opts = optiset('solver','mkltrnls','display','iter','solverOpts',mkltrnlsset('tolTR',1e-9));
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'data',xdata,ydata,'ndec',3,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt)

%% NLS3 (Modified NLP - HS76)
clc
%Function
fun = @(x) [x(1);
            sqrt(0.5)*x(2);
            x(3);
            sqrt(0.5)*x(4);];
%Fitting Data
ydata = [0.0, 0.0, 0.0, 0.0];
%Constraints
lb = [0.0, 0.0, 0.0, 0.0];
ub = [inf, inf, inf, inf];
A = -[-1.0, -2.0, -1.0, -1.0;
     -3.0, -1.0, -2.0, 1.0];
b = -[-5.0, -0.4]';
Aeq = [0.0, 1.0, 4.0, 0.0];
beq = 1.5;
%Setup Options
opts = optiset('solver','levmar','display','iter');
%Build & Solve
x0 = [0.5, 0.5, 0.5, 0.5];
Opt = opti('fun',fun,'ydata',ydata,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)

%% NLS Example Prob
clc
%Get Problem
prob = nls_prob(19);
opts = optiset('display','iter');
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
opts = optiset('solver','nlopt','display','iter');
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

%% MINLP1
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
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'ineq',A,b,'eq',Aeq,beq)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% MINLP2
clc
%Objective
obj = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);
% Constraints
A = [-1 2.5; 1 2.5]; 
b = [1;-15];
% Solve
Opt = opti('obj',obj,'ndec',2,'ivars',[1 2],'ineq',A,b)
x0 = [0;0];
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% MINLP3
clc
%Objective
obj = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
% Constraints
nlcon = @(x) -x(1)^2 + x(2)-0.5;
nlrhs = 0;
nle = 1;      
% Solve
Opt = opti('obj',obj,'ndec',2,'nlmix',nlcon,nlrhs,nle,'int',[1 2],'options',optiset('solver','scip','display','iter'))
x0 = [4.9;0.1];
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[])

%% BONMIN Example MINLP
clc
%Objective
obj = @(x) -x(1) - x(2) - x(3);
%Constraints
nlcon = @(x) [ (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.);
                x(1) - x(2);
                x(1) + x(3) + x(4)];
nlrhs = [1/4;0;2];
nle = [-1;-1;-1];
ub = [1;Inf;Inf;5];
lb = [0;0;0;0];
% Solve
bopts = bonminset('algorithm','B-OA','milp_solver','Cbc_d');
opts = optiset('solver','bonmin','display','iter','solverOpts',bopts);
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'int','BCCI','options',opts)
x0 = [0;0;0;0];
[x,fval,ef,info] = solve(Opt,x0)

%% MATLAB GA Example MINLP
clc
%Objective
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
%Constraints
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];
%Setup Options
opts = optiset('solver','nomad','display','iter');
% Solve
Opt = opti('obj',fun,'bounds',lb,ub,'int','IC','options',opts)
x0 = [16;0];
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,1)

%% Large Sparse System of Linear Equations
n = 3000;
A = sprandn(n,n,0.2);
b = randn(n,1);

% Build OPTI Problem
Opt = opti('sle',A,b,'options',optiset('solver','mumps'))
Opt2 = opti('sle',A,b,'options',optiset('solver','matlab'))
 
% Solve
tic
x = solve(Opt);
toc

tic
x2 = solve(Opt2);
toc

norm(x-x2)

%% Large Complex Sparse System of Linear Equations
n = 1500;
A = sprandn(n,n,0.2) + sparse(1:n,1:n,1i*ones(n,1));
b = randn(n,1);

% Solve Directly (OPTI doesn't support complex numbers)
tic
x = opti_zmumps(A,b);
toc

tic
x2 = A\b;
toc

norm(x-x2)


