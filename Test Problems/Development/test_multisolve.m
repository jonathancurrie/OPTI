%% Nonconvex QP
clc
%Objective & Constraints
H = sparse([0 -1; -1 0]);
f = [0;0];
%Constraints
lb = [-0.5;-0.5];
ub = [1;1];

%Build & Solve
Opt = opti('qp',H,f,'bounds',lb,ub,'options',optiset('solver','ooqp','display','final','warnings','critical'))
[x,fval,exitflag,info] = multisolve(Opt)
%Plot
figure(1)
plot(Opt)
figure(2)
multiplot(Opt)

%% Nonconvex QP w Inequalities and Maximization
clc
%Objective & Constraints
H = sparse([0 -1; -1 0]);
f = [0;0];
%Constraints
A = [1 1; 1 -1]; b = [1.5;1.2];
lb = [-0.5;-0.5];
ub = [1;1];

%Build & Solve
Opt = opti('qp',H,f,'ineq',A,b,'bounds',lb,ub,'sense',-1,'options',optiset('solver','ooqp','display','iter','warnings','all'))
[x,fval,exitflag,info] = multisolve(Opt)
%Plot
figure(1)
plot(Opt)
figure(2)
multiplot(Opt)


%% Nonconvex QP2
clc
%Objective & Constraints
H = sparse([0 0 -1;0 -1 0;-1 0 0]);
f = [0;0;0];
%Constraints
lb = [0;0;0]; ub = [4;4;4];

%Build & Solve
Opt = opti('qp',H,f,'bounds',lb,ub,'options',optiset('solver','auto','display','final'))
[x,fval,exitflag,info] = multisolve(Opt)
multiplot(Opt)

%% NLS1
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Setup Options
opts = optiset('solver','mkltrnls','display','final');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'bounds',[400;-1],[600;1],'options',opts)
[x,fval,exitflag,info] = multisolve(Opt)
%Plot
figure(1)
plot(Opt)
figure(2)
multiplot(Opt,1)


%% NLS2
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Setup Options
opts = optiset('solver','mkltrnls','display','final');
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'ydata',ydata,'bounds',[0;0;0],[10;5;5],'options',opts,'x0',x0)
[xs,fvals] = multisolve(Opt,x0)

plot(Opt)


%% NLP1
clc
% Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Linear Constraints
A = [-1 1];  b = -1;
Aeq = [1.1 1]; beq = 5; 
lb = [0;0]; ub = [4;4];
% Initial Guess
x0 = [2;2];

% Build OPTI Object
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ineq',A,b,'eq',Aeq,beq,'options',optiset('display','final'))
% Solve with Violation Penalty of 1e2
[x,fval,ef,info] = multisolve(Opt,x0)
%Plot
multiplot(Opt,1)

%% NLP2
clc
% Objective
 fun = @(x) norm([x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
                  x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))]);

% Bounds
 lb = [-4;-4]; 
 ub = [4;4];

% Initial Solution
 x0 = [-4;-4];

% Build optiprob structure (intermediate structure)
 Opt = opti('fun',fun,'bounds',lb,ub,'x0',x0,'options',optiset('solver','filterSD','display','iter'));
 
 [x,fval,ef,info] = multisolve(Opt,x0)
%Plot
figure(1)
plot(Opt,[],1)
figure(2)
multiplot(Opt,1)

%% NLP3
 % Objective
 fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));

 A = [1 1]; b = -1.1;
 
% Bounds
 lb = [-5;-5];
 ub = [5;5];

%Initial Guess
 x0 = [0;0];
 
 % Build optiprob structure (intermediate structure)
 Opt = opti('fun',fun,'ineq',A,b,'bounds',lb,ub,'x0',x0,'options',optiset('solver','filterSD','display','iter'));
 
 [x,fval,ef,info] = multisolve(Opt,x0)
%Plot
figure(1)
plot(Opt,[-2 0 -2 0])
figure(2)
multiplot(Opt)