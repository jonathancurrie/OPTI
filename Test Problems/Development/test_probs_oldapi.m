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
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'bounds',[0;0],[10;10]);
opts = optiset('solver','auto'); %automatically choose the best solver available
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)  
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'eq',[1 1 1],40,'bounds',[0 0 0]',[40 inf inf]');
opts = optiset('solver','qsopt');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% LP3
clc
%Objective & Constraints
f = [8,1]';
A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
b = [-4,2,21,39,3,0,0]';
%Build & Options
prob = optiprob('grad',f,'ineq',A,b);
opts = optiset('solver','clp');
%Solve
x0 = [9,6]';
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,10)

%% SLE
clc
clear
%Problem
n = 800;
A = sprandn(n,n,n*n/10);
b = 1:n; b = b';
%Build & Options
prob = optiprob('sle',A,b);
opts = optiset('solver','mumps');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt);
info

%% MILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];  
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II');
opts = optiset('solver','glpk');
%Solve
Opt = opti(prob,opts)
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
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',[0 0 0 2]',[40 inf inf 3]','int','CCCI');
opts = optiset('solver','cbc');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% MILP3 Infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = [4,1; 1.5,1; -2, 1; -0.2, 1]; 
b = [5.3;4.5;-2.5; 1.5];  
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II');
opts = optiset('solver','lp_solve');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP4
clc
%Objective & Constraints
f = -[-1, 2]';
A = [2, 1;-4, 4];
b = [5, 5]';
e = -[1, 1];
%Build & Options
prob = optiprob('grad',f,'mix',A,b,e,'int','II');
opts = optiset('solver','auto');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP5
clc
%Objective & Constraints
f = [3, -7, -12]';
A = [-3, 6, 8;6, -3, 7;-6, 3, 3];
b = [12, 8, 5]';
e = [-1, -1, -1];
%Build & Options
prob = optiprob('grad',f,'mix',A,b,e,'int','III');
opts = optiset('solver','auto');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% BILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'int','BB');
opts = optiset('solver','glpk');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,3)

%% BILP2
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
b = [9; 1; 0; 0];
%Build & Options
prob = optiprob('grad',f,'ineq',A,b,'int','BBBB');
opts = optiset('solver','glpk');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% QP1
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Build & Options
prob = optiprob('hess',H,'grad',f,'ineq',A,b);
opts = optiset('solver','ooqp');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% QP2
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
%Build & Options
prob = optiprob('qp',H,f,'ineq',A,b,'lb',[0;0]);
opts = optiset('solver','clp');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
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
%Build & Options
prob = optiprob('qp',H,f,'eq',Aeq,beq,'ineq',A,b,'bounds',[0;0],[10;10]);
opts = optiset('solver','clp');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% QP4 unconstrained
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
prob = optiprob('hess',H,'grad',f);
%Solve
Opt = opti(prob)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% QCQP1
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
prob = optiprob('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub);
%Solve
Opt = opti(prob)
[x,fval,exitflag,info] = solve(Opt)

%% QCQP2
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
prob = optiprob('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub);
%Solve
Opt = opti(prob)
[x,fval,exitflag,info] = solve(Opt)
% Plot
plot(Opt)

%% MIQP1
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
b = [3; 3.5];
ivars = 1;
%Build & Options
prob = optiprob('hess',H,'grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int',ivars);
opts = optiset('solver','auto');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,4)

%% MIQP2
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Build & Options
prob = optiprob('hess',H,'grad',f,'ineq',A,b,'int','CIC');
opts = optiset('solver','auto');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% MIQCQP1
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
%Build & Options
prob = optiprob('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'int',int);
opts = optiset('solver','auto');
% Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)

%% SNLE1
clc
%Equations
f = @(x) [2*x(1) - x(2) - exp(-x(1));
          -x(1) + 2*x(2) - exp(-x(2))];
%Build & Options
x0 = [-5;5];
prob = optiprob('fun',f,'x0',x0);
Opt = opti(prob)
%Solve
[x,fval,ef,stat] = solve(Opt)

%% NLS1
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Build & Options
prob = optiprob('fun',fun,'data',xdata,ydata,'ndec',2);
opts = optiset('solver','lmder','display','iter');
%Solve
x0 = [100; -1]; % Starting guess
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,0.05,1)

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
%Build & Options
prob = optiprob('fun',fun,'ydata',ydata,'ndec',3);
opts = optiset('solver','nl2sol');
%Solve
x0=[1.0; 0.0; 0.0];
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt,x0)

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
%Build & Options
prob = optiprob('fun',fun,'ydata',ydata,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub);
opts = optiset('solver','levmar');
%Solve
x0 = [0.5, 0.5, 0.5, 0.5];
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt,x0)

%% UNO Rosenbrock
clc
%Objective & Gradient
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;

%Build & Options
prob = optiprob('obj',obj,'ndec',2);
opts = optiset('solver','nlopt');
%Solve
x0 = [0 0]';
Opt = opti(prob,opts)
[x,fval,exitflag,info]= solve(Opt,x0)
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
%Build & Options
prob = optiprob('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub); %ndec not required if we can determine ndec from bounds
opts = optiset('solver','ipopt','warnings','on','display','iter');
%Solve
x0 = [1 5 5 1]';
Opt = opti(prob,opts)
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
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];
%Hessian (must not be tril for matlab!) & Hessian Structure       
hess = @(x) sparse(  [ 1200*x(1)^2-400*x(2)+2  -400*x(1)       0                          0
                       -400*x(1)               220.2           0                          19.8
                        0                      0               1080*x(3)^2- 360*x(4) + 2  -360*x(3)
                        0                      19.8            -360*x(3)                   200.2 ]);  
Hstr = @() sparse([ 1  0  0  0 
                    1  1  0  0
                    0  0  1  0
                    0  1  1  1 ]);                               
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Build & Options
prob = optiprob('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub);
opts = optiset('solver','ipopt','maxfeval',1e7,'maxiter',1e6);
%Solve
x0 = [-3  -1  -3  -1]';
Opt = opti(prob,opts)
[x,fval,exitflag,info]= solve(Opt,x0)

%% NLP3 Hock & Schittkowski #51
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
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
%Build & Options
prob = optiprob('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr);
opts = optiset('solver','nlopt');
%Solve
x0 = [ 2.5 0.5 2 -1 0.5 ];
Opt = opti(prob,opts)
[x,fval,exitflag,info]= solve(Opt,x0)

%% NLP4 Hock & Schittkowski #104 (not working)
clc
%Objective
obj = @(x) 0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
%Constraints
nlcon = @(x) [ 1 - 0.588*x(5)*x(7) - 0.1*x(1);
               1 - 0.588*x(6)*x(8) - 0.1*x(1) - 0.1*x(2);
               1 - 4*x(3)*x(5)^(-1) - 2*x(3)^(-0.71)*x(5)^(-1) - 0.588*x(3)^(-1.3)*x(7);
               1 - 4*x(4)*x(6)^(-1) - 2*x(4)^(-0.71)*x(6)^(-1) - 0.588*x(4)^(-1.3)*x(8);
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2)];
nlrhs = [ 0 0 0 0 1 4.2]';
nle = [1 1 1 1 1 -1]';
lb = 0.1*ones(8,1);
ub = 10*ones(8,1);
%Build & Options
prob = optiprob('obj',obj,'nlmix',nlcon,nlrhs,nle,'nljac',symJac(nlcon),'bounds',lb,ub);
opts = optiset('solver','ipopt','maxiter',50000,'display','iter');
%Solve
x0 = [6,3,.4,.2,6,6,1,.5];
Opt = opti(prob,opts)
[x,fval,exitflag,info]= solve(Opt,x0)

%% NLP HS Interface
clc
%Get HS problem
no = 15;
[prob,sol,fmin] = nlp_HS(no);
%Build & Options
opts = optiset('solver','ipopt');
%Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info]= solve(Opt)

acc = norm(fmin-fval)

plot(Opt,[],1) %problems 1-24 can be plotted

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
%Build & Options
prob = optiprob('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'ineq',A,b,'eq',Aeq,beq);
% Solve
x0 = [2;2];
Opt = opti(prob)
[x,fval,ef,info] = solve(Opt,x0)

plot(Opt,[],1)

%% MINLP2
clc
%Objective
obj = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);
% Constraints
A = [-1 2.5; 1 2.5]; 
b = [1;-15];
%Build & Options
prob = optiprob('obj',obj,'ndec',2,'ivars',[1 2],'ineq',A,b);
% Solve
x0 = [0;0];
Opt = opti(prob)
[x,fval,ef,info] = solve(Opt,x0)

plot(Opt)

%% MINLP3
clc
%Objective
obj = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
% Constraints
nlcon = @(x) -x(1)^2 + x(2)-0.5;
nlrhs = 0;
nle = 1;      
%Build & Options
prob = optiprob('obj',obj,'ndec',2,'nlmix',nlcon,nlrhs,nle,'int',[1 2]);
% Solve
x0 = [4.9;0.1];
Opt = opti(prob)
[x,fval,ef,info] = solve(Opt,x0)

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
%Build & Options
prob = optiprob('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'int','BCCI');
% Solve
x0 = [0;0;0;0];
Opt = opti(prob)
[x,fval,ef,info] = solve(Opt,x0)

%% MATLAB GA Example MINLP
clc
%Objective
obj = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
%Constraints
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];
%Build & Options
prob = optiprob('obj',obj,'bounds',lb,ub,'int','IC');
opts = optiset('solver','nomad','display','iter');
% Solve
x0 = [16;0];
Opt = opti(prob,opts)
[x,fval,ef,info] = solve(Opt,x0)

plot(Opt,1)

%% Gradients
clear
clc

syms x1 x2

%OPTI Sym
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1
jac = symJac(fun)

%Symbolic
sfun = sin(x1+x2) + (x1-x2)^2 - 1.5*x1 + 2.5*x2 + 1
jac1 = jacobian(sfun)

%OPTI Sym
fun = @(x) [x(1) + x(2) - 1;      
            x(1)^2 + x(2)^2 - 1];
jac = symJac(fun)

%Symbolic
sfun = 0.01*(x1-1)^2 + (x2-x1^2)^2;
jac1 = jacobian(sfun)











