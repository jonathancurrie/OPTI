%% Dense vs Sparse SNLE
clc
fun = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('fun',fun,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)


%% Sparse above
fun = @(x) -1;
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
cl = zeros(4,1);
cu = zeros(4,1);
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('fun',fun,'nl',nleq,cl,cu,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [no grad]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad DENSE]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
[nljac,nljacstr] = symJac(nleq);        
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'nljac',nljac,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad SPARSE]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
[nljac,nljacstr] = symJac(nleq); nljac = @(x) sparse(nljac(x));       

x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'nljac',nljac,'nljacstr',nljacstr,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad SPARSE alt nl format]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
[nljac,nljacstr] = symJac(nleq); nljac = @(x) sparse(nljac(x));       
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'nljac',nljac,'nljacstr',nljacstr,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad SPARSE alt nl format w lin]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
[nljac,nljacstr] = symJac(nleq); nljac = @(x) sparse(nljac(x));       
A = [0 1 0 1]; b = 2;
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'ineq',A,b,'nljac',nljac,'nljacstr',nljacstr,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad SPARSE alt nl format w lin + int]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];
[nljac,nljacstr] = symJac(nleq); nljac = @(x) sparse(nljac(x));       
A = [0 1 0 1]; b = 2;
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nleq',nleq,'ivars',2,'lin',A,-Inf,b,'nljac',nljac,'nljacstr',nljacstr,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% As above but new construct [w grad SPARSE alt nl format w ineq]
clc
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4)-x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))
            x(3)];
[nljac,nljacstr] = symJac(nleq); nljac = @(x) sparse(nljac(x));       
nlrhs = zeros(5,1);
nle = [zeros(4,1);-1];
A = [0 1 0 1]; b = 2;
        
x0 = [-30 -10 -30 -10]';        

Opt = opti('nlmix',nleq,nlrhs,nle,'ineq',A,b,'nljac',nljac,'nljacstr',nljacstr,'x0',x0,'options',optiset('solver','auto'))

[x,f,e,i] = solve(Opt)

%% Wiki Ex 1
clc
% System of Nonlinear Equations
 nleq = @(x) [ 2*x(1) - x(2) - exp(-x(1));
             -x(1) + 2*x(2) - exp(-x(2))];

% Starting Guess
 x0 = [-5;5];

% Create OPTI Object
 Opt = opti('nleq',nleq,'x0',x0)

% Solve the SNLE problem
[x,fval,exitflag,info] = solve(Opt)

%% Online Example
clc
% System of Nonlinear Equations
 nleq = @(x) [10*(x(2) - x(1)^2)
             sqrt(90)*(x(4) - x(3)^2)
             sqrt(10)*(x(2) + x(4) - 2)
             (1/sqrt(10))*(x(2) - x(4))];

% Nonlinear Equations Jacobian
 nlJac = @(x) sparse([-20*x(1),10,0,0
                      0,0,-6*10^(1/2)*x(3),3*10^(1/2)
                      0,10^(1/2),0,10^(1/2)
                      0,10^(1/2)/10,0,-10^(1/2)/10]);

% Jacobian Sparsity Pattern
 nlJacstr = @() sparse([1 1 0 0
                        0 0 1 1
                        0 1 0 1
                        0 1 0 1]);

% Starting Guess
 x0 = [-30;-10;-30;-10];

% Sparse SNLE OPTI Problem
 Opt = opti('nleq',nleq,'nlJac',nlJac,'nlJacstr',nlJacstr,'x0',x0)

  %Solve
[x,f,e,i] = solve(Opt)

%% SCNLE Example
clc

% Objective (Nonlinear Equations) Function
 fun = @(x) [ 2*x(1) - x(2) - exp(-x(1));
             -x(1) + 2*x(2) - exp(-x(2))];

% Bounds
 lb = [0.6;0];
 ub = [1;1];

% Starting Guess
 x0 = [-5;5];

% Create OPTI Object
 Opt = opti('nleq',fun,'bounds',lb,ub,'x0',x0,'options',optiset('display','iter'))

% Solve the SCNLE problem
[x,fval,exitflag,info] = solve(Opt)


%% Problem SCNLE
clc

A = randn(91);
B12 = randn(91,145);
B1 = randn(145,91); B2 = randn(145,91);
Constant = randn(91,1);
A_ieq = randn(12,91); b_ieq = randn(12,1);
lb = zeros(91,1); ub = 100*ones(91,1);
z0 = lb;

nleq = @(x) A*x + B12*((B1*x).*(B2*x)) + Constant; 

Opt=opti('nleq',nleq,'ineq',A_ieq,b_ieq,'bounds',lb,ub,'x0',z0,'options',optiset('display','iter'))

solve(Opt)
