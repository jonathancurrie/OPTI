%% LP1
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
%Build Object
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10]);
% [x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','qsopt');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1 1],40,'bounds',[0 0 0]',[40 inf inf]','options',opts)
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

%% LP3
clc
%Objective & Constraints
f = [8,1]';
A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
b = [-4,2,21,39,3,0,0]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
x0 = [9,6]';
Opt = opti('grad',f,'ineq',A,b,'options',opts,'x0',x0)
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

%% QP2
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('qp',H,f,'ineq',A,b,'lb',[0;0],'options',opts)
[x,fval,exitflag,info] = solve(Opt)
% write(Opt,'test.gms')

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
write(Opt,'test.gms')

%% QP4 unconstrained
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
%Build & Solve
Opt = opti('hess',H,'grad',f)
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

%% Nonconvex QP
clc
%Objective & Constraints
H = sparse([0 -1; -1 0]);
f = [0;0];
%Constraints
lb = [-0.5;-0.5];
ub = [1;1];

%Build & Solve
Opt = opti('qp',H,f,'bounds',lb,ub,'options',optiset('solver','scip','display','iter'))
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

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
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('display','iter','solver','scip'))
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

%% Positive Semidefinite QCQP
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
write(Opt,'test.gms')

%% Indefinite QCQP
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
write(Opt,'test.gms')

%% MIQP1 [-11.5]
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
b = [3; 3.5];
ivars = 1;
%Setup Options
opts = optiset('display','iter','solver','auto');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int',ivars,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

%% MIQP2 [-2.75]
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','cplex');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'int','CIC','options',opts)
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'test.gms')

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
write(Opt,'test.gms')

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
write(Opt,'test.gms')

%% QCQP Equals
clc
% Objective
 H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
 f = -[2 2]';                
% Linear Constraints
 A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
 b = [2;5];    
 lb = [0;0];                 %Bounds on x (lb <= x)
% Quadratic Constraints
 Q = {[1 0; 0 1]             %Quadratic Constraints (qrl <= x'Qx + l'x <= qru)
      [1 0; 0 1]};
 l = {[0;-2]; [-2;2]};
 qrl = {3; 1};               %QC1 is double sided, QC2 is an equality
 qru = {5; 1};
% Create OPTI Object
 Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qcrow',Q,l,qrl,qru,'options',optiset('solver','ipopt'))
% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)
write(Opt,'ncQCQP.gms')
plot(Opt)

%% UNO 1
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','nlopt','display','off');
%Build & Solve
x0 = [0 0]';

Opt = opti('obj',obj,'ndec',2,'options',opts,'x0',x0)
[x,fval,exitflag,info]= solve(Opt)
write(Opt,'test.gms')

%% NLP1 Hock & Schittkowski #71
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
lb = ones(4,1); lb(3) = -Inf;
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];         
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','all','display','iter');
%Build & Solve
x0 = [1 5 5 1]';
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts,'x0',x0)
[x,fval,exitflag,info]= solve(Opt)
write(Opt,'test.gms')

%% NLP2 Hock & Schittkowski #38
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);                              
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Setup Options
opts = optiset('solver','ipopt');
%Build & Solve
x0 = [-3  -1  -3  -1]';
Opt = opti('obj',obj,'bounds',lb,ub,'options',opts,'x0',x0)
[x,fval,exitflag,info]= solve(Opt)
write(Opt,'test.gms')

%% NLP3 Hock & Schittkowski #51
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
%Nonlinear Constraints & Jacobian Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nlrhs = [4 0 0]';
nle = [0 0 0]';
%Setup Options
opts = optiset('solver','nlopt');
%Build & Solve
x0 = [ 2.5 0.5 2 -1 0.5 ];
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'options',opts,'x0',x0)
[x,fval,exitflag,info]= solve(Opt)
write(Opt,'test.gms')

%% NLP4 Hock & Schittkowski #104 (not working)
clc
%Objective
obj = @(x) 10 + 0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) - x(1) - x(2);
%Constraints
nlcon = @(x) [ -0.588*x(5)*x(7) - 0.1*x(1);
               -0.588*x(6)*x(8) - 0.1*x(1) - 0.1*x(2);
               -4*x(3)/x(5) - 2/(x(3)^(0.71)*x(5)) - 0.588*x(7)/(x(3)^(1.3));
               -4*x(4)/x(6) - 2/(x(4)^(0.71)*x(6)) - 0.588*x(8)/(x(4)^(1.3));
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) - x(1) - x(2);
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) - x(1) - x(2)];
nlrhs = [ -1 -1 -1 -1 -9.9 -5.8]';
nle = [1 1 1 1 1 -1]';
lb = 0.1*ones(8,1);
ub = 10.1*ones(8,1);
%Setup Options
opts = optiset('solver','scip','display','iter');
%Build & Solve
x0 = [6,3,.4,.2,6,6,1,.5];
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts,'x0',x0)
[x,fval,exitflag,info]= solve(Opt)
write(Opt,'test.gms')

%% MINLP1
clc
% Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Linear Constraints
A = [-1 1]; b = -1;
Aeq = [1 1]; beq = 5; 
lb = [0;0]; ub = [4;4];
%Integer Constraints
xtype = 'IC';
Opt = opti('obj',obj,'bounds',lb,ub,'xtype',xtype,'ineq',A,b,'eq',Aeq,beq)
write(Opt,'minlp.gms')

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
write(Opt,'test.gms')

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
bopts = bonminset('algorithm','B-OA','milp_solver','Cplex');
opts = optiset('display','iter','solverOpts',bopts);
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'int','BCCI','options',opts)
x0 = [0;0;0;0];
[x,fval,ef,info] = solve(Opt,x0)
write(Opt,'test.gms')
