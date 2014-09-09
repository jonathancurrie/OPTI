%% LP1
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
lb = [0;0]; ub = [10;10];
%Setup Options
opts = optiset('solver','clp'); 
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt)  
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% LP1 ROW
rl = -Inf(3,1);
ru = b;
Opt = opti('grad',f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt)  
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% Raw Test
opts.debug = 1;
opts.display = 2;
opts.maxiter = 5;
% rl = zeros(3,1);
n = length(f);
H = spalloc(n,n,0);
Aeq = spalloc(0,n,0);
beq = zeros(0,1);

[x,fv,ef,iter,l] = ooqp([],f,sparse(A),rl,ru,[],[],lb,ub,opts)

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
Aeq = [1 1 1]; beq = 40;
lb = [0 0 0]'; ub = [40 inf inf]';
%Setup Options
opts = optiset('solver','ooqp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
    
%% LP2 ROW
Ar = [A;Aeq];
rl = [-Inf(2,1); beq];
ru = [b;beq];
Opt = opti('grad',f,'lin',Ar,rl,ru,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% LP2 ROW + EQ (OOQP FORMAT) (NOT SUPPORTED CURRENTLY)
% rl = -Inf(2,1);
% ru = b;
% prob = optiprob('grad',f,'lin',A,rl,ru,'eq',Aeq,beq,'bounds',lb,ub);
% Opt = opti(prob,optiset('solver','clp'))
% [x,fval,exitflag,info] = solve(Opt)  
% %Check Solution
% [ok,msg] = checkSol(Opt)

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
Opt = opti('grad',f,'ineq',A,b,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)
plot(Opt,10)

%% LP3 ROW
rl = -Inf(7,1);
ru = b;
Opt = opti('grad',f,'lin',A,rl,ru,'options',opts)
[x,fval,exitflag,info] = solve(Opt)  
%Check Solution
[ok,msg] = checkSol(Opt)


%% MILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];  
lb = [0;0]; ub=[10;10];
%Setup Options
opts = optiset('solver','glpk');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',lb,ub,'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)


%% MILP1 ROW
rl = -Inf(3,1);
ru = b;
Opt = opti('grad',f,'lin',A,rl,ru,'bounds',lb,ub,'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)  
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)


%% MILP2
clc
%Objective & Constraints
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
lb = [0 0 0 2]';
ub = [40 inf inf 3]';
int = 'CCCI';
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'int',int,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% MILP2 ROW
Ar = [A;Aeq];
rl = [-Inf(2,1); beq];
ru = [b;beq];
Opt = opti('grad',f,'lin',Ar,rl,ru,'bounds',lb,ub,'int',int,'options',opts)
[x,fval,exitflag,info] = solve(Opt)  
%Check Solution
[ok,msg] = checkSol(Opt)


%% BILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BB','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,3)

%% BILP1 ROW
rl = -Inf(4,1);
ru = b;
Opt = opti('grad',f,'lin',A,rl,ru,'int','BB','options',opts)
[x,fval,exitflag,info] = solve(Opt)  
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% BILP2
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
b = [9; 1; 0; 0];
%Setup Options
opts = optiset('solver','cbc');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BBBB','options',opts)
[x,fval,exitflag,info] = solve(Opt)

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
opts = optiset('solver','matlab','warnings','on','display','iter');
%Build & Solve
x0 = [1 5 5 1]';
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
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
               x(3) + x(4) - 2*x(5)];
%                x(2) - x(5)];
nljac = @(x) sparse([ 1  3  0  0  0;
	                  0  0  1  1 -2]);
% 	                  0  1  0  0 -1]);
nljacstr = @() sparse([1 1 0 0 0;
                       0 0 1 1 1]);
%                        0 1 0 0 1]);
                   
cl = [4;0];
cu = [4;0];
nlrhs = [4 0 1 15]';
nle = [0 0 1 -1]';

A = [0 1 0 0 -1];
rl = 1;
ru = 15;

%Setup Options
opts = optiset('solver','ipopt','display','iter','maxiter',1e4);
%Build & Solve
x0 = [ 2.5 0.5 2 -1 2.5 ];
Opt = opti('obj',obj,'grad',grad,'lin',A,rl,ru,'nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'options',opts)
[x,fval,exitflag,info]= solve(Opt,x0)

info.Lambda
info.Lambda.ineqlin
info.Lambda.ineqnonlin
info.Lambda.eqnonlin

%% Test Conversion Function

Opt.prob.nlcon(x0)
full(Opt.prob.nljac(x0))
full(Opt.prob.nljacstr())



%% SINGLE BOUNDS NO JAC
clc
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5)];

cl = [-Inf 1 2];
cu = [4 1 Inf];

x0 = [ 2.5 0.5 2 -1 1.5 ];

prob = optiprob('nl',nlcon,cl,cu,'x0',x0)
prob = nrow2mix(prob)

prob.nlrhs
prob.nle
prob.nlcon(x0)

%% SINGLE BOUNDS W JAC
clc
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5)];

nljac = @(x) sparse([ 1 3 0 0 0;
                      0 0 1 1 -2;
                      0 1 0 0 -1]);

cl = [-Inf 1 2];
cu = [4 1 Inf];

x0 = [ 2.5 0.5 2 -1 1.5 ];

prob = optiprob('nl',nlcon,cl,cu,'nljac',nljac,'x0',x0)
prob = nrow2mix(prob)

prob.nlrhs
prob.nle
prob.nlcon(x0)

%% DUAL BOUNDS W JAC
clc
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5);
               x(1) + 3*x(2);
               x(3) + x(4) - 5*x(5)];

nljac = @(x) sparse([ 1 3 0 0 0;
                      0 0 1 1 -2;
                      0 1 0 0 -1;
                      1 3 0 0 0;
                      0 0 1 1 -5;]);

cl = [-Inf 1 -Inf 2 5];
cu = [4 1 12 16 Inf];

x0 = [ 2.5 0.5 2 -1 1.5 ];

prob = optiprob('nl',nlcon,cl,cu,'nljac',nljac,'x0',x0)
prob = nrow2mix(prob)

prob.nlrhs
prob.nle
prob.nlcon(x0)
full(prob.nljac(x0))

%% DUAL BOUNDS W JAC & STR
clc
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5);
               x(1) + 3*x(2);
               x(3) + x(4) - 5*x(5)];

nljac = @(x) sparse([ 1 3 0 0 0;
                      0 0 1 1 -2;
                      0 1 0 0 -1;
                      1 3 0 0 0;
                      0 0 1 1 -5;]);

nljacstr = @() sparse([ 1 1 0 0 0; 
                        0 0 1 1 1;
                        0 1 0 0 1;
                        1 1 0 0 0;
                        0 0 1 1 1]);

cl = [-Inf 1 -Inf 2 5];
cu = [4 1 12 16 Inf];

x0 = [ 2.5 0.5 2 -1 1.5 ];

prob = optiprob('nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'x0',x0)
prob = nrow2mix(prob)

prob.nlrhs
prob.nle
prob.nlcon(x0)
full(prob.nljac(x0))
full(prob.nljacstr())


%% Ali Constraint Linearity Index Problem
clc
obj=@(x) (x)^2;
A=1; ru=1; rl=-1; xtype='I'; x0=1;
opt=opti('obj',obj,'xtype',xtype,'x0',x0,'lin',A,rl,ru)
opt.prob.sizes
solve(opt,x0)

%% Identical form of above
clc
obj=@(x) (x)^2;
A=[1;-1];  b=[1;1];  xtype='I'; x0=1;
opt=opti('obj',obj,'xtype',xtype,'x0',x0,'ineq',A,b); 
opt.prob.sizes
solve(opt,x0)

%% Nonlinear Constraints Double Sided Row, Linear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
A = [-1 1; 1 1];
rl = [-3;5];
ru = [-1;5];
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'lin',A,rl,ru)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% Nonlinear Constraints Double Sided Ineq, Linear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
A = [-1 1; 1 -1];
b = [-1;3];
Aeq = [1 1];
beq = 5;
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'ineq',A,b,'eq',Aeq,beq)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% Nonlinear Constraints Double Sided Row, NonLinear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
nlcon = @(x) [-x(1) + x(2); x(1) + x(2)];
cl = [-3;5];
cu = [-1;5];
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'nl',nlcon,cl,cu)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% Nonlinear Constraints Double Sided Ineq, Nonlinear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
nlcon = @(x) [-x(1) + x(2); x(1) - x(2); x(1) + x(2)];
nlrhs = [-1;3;5];
nle = [-1;-1;0];
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'nlmix',nlcon,nlrhs,nle)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% Nonlinear Constraints Double Sided Ineq, Nonlinear + Linear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
A = [-1 1; 1 -1];
b = [0;2];
nlcon = @(x) [-x(1) + x(2); x(1) - x(2); x(1) + x(2)];
nlrhs = [-1;3;5];
nle = [-1;-1;0];
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'nlmix',nlcon,nlrhs,nle,'ineq',A,b)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)

%% Nonlinear Constraints Double Sided Row, NonLinear + Linear
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Constraints
A = [-1 1];
rl = -2; ru = 0;
nlcon = @(x) [-x(1) + x(2); x(1) + x(2)];
cl = [-3;5];
cu = [-1;5];
lb = [0;0]; ub = [4;4];
% Solve
x0 = [2;2];
Opt = opti('obj',obj,'ndec',2,'bounds',lb,ub,'ivars',1,'nl',nlcon,cl,cu,'lin',A,rl,ru)
[x,fval,ef,info] = solve(Opt,x0)
%Plot
plot(Opt,[],1)
