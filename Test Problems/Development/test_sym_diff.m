
%% Gradients
clear
clc

syms x1 x2
x = [0.1;0.1];

%OPTI Sym
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1
[jac,jacp] = symJac(fun)
[hess,hessp] = symHess(fun)

%Symbolic
sfun = sin(x1+x2) + (x1-x2)^2 - 1.5*x1 + 2.5*x2 + 1
jac1 = jacobian(sfun)
hess1 = hessian(sfun)

%Eval
[jac(x);double(subs(jac1,{'x1','x2'},{x(1) x(2)}))]
full(jacp())
[hess(x);double(subs(hess1,{'x1','x2'},{x(1) x(2)}))]
full(hessp())

%OPTI Sym
fun = @(x) [x(1) + x(2) - 1;      
            x(1)^2 + x(2)^2 - 1];
jac = symJac(fun)

%Symbolic
fun = [x1 + x2 - 1;      
       x1^2 + x2^2 - 1];
jac1 = jacobian(fun)

%Eval
[jac(x);double(subs(jac1,{'x1','x2'},{x(1) x(2)}))]

%OPTI Sym
fun = @(x) 0.01*(x(1)-1)^2 + (x(2)-x(1)^2)^2;
[jac,jacp] = symJac(fun)
[hess,hessp] = symHess(fun)

%Symbolic
sfun = 0.01*(x1-1)^2 + (x2-x1^2)^2;
jac1 = jacobian(sfun)
hess1 = hessian(sfun)

%Eval
[jac(x);double(subs(jac1,{'x1','x2'},{x(1) x(2)}))]
full(jacp())
[hess(x);double(subs(hess1,{'x1','x2'},{x(1) x(2)}))]
full(hessp())


%% Symbolic Hessian of Lagrangian
clc
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);
% Nonlinear Constraints 
 nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
                x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];

jac = symJac(fun)          
hessObj = symHess(fun)

[hessLag,pattern] = symHessLag(fun,nlcon,4,true)


% Hessian
 H = @(x,sigma,lambda) sparse(sigma*[ 2*x(4)             0      0   0;
                                      x(4)               0      0   0;
                                      x(4)               0      0   0;
                                      2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ] + ...
                      lambda(1)*[   0          0         0         0;
                                 x(3)*x(4)     0         0         0;
                                 x(2)*x(4) x(1)*x(4)     0         0;
                                 x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ] + ...
                      lambda(2)*diag([2 2 2 2]));
                  
                  
x = 0.1:0.1:0.4;
sigma = 0.3;
lambda = 1:2;

hessLag(x,sigma,lambda)
full(H(x,sigma,lambda))
full(pattern())


%% Application to an NLP, first no derivs
clc
% Objective
 fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Nonlinear Constraints 
 nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
                x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
 cl = [25;40];
 cu = [Inf;40]; 

% Bounds (lb <= x <= ub)
 lb = ones(4,1);
 ub = 5*ones(4,1);         

% Initial Guess
 x0 = [1 5 5 1]';
 
opts = optiset('display','iter','solver','ipopt','derivCheck','on','solverOpts',ipoptset('derivative_test','first-order'));
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'x0',x0,'options',opts)

[x,f,e,i] = solve(Opt)


%% Application to an NLP, now with symJac, symHessLag
clc
% Objective
 fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

grad = symJac(fun); 
 
% Nonlinear Constraints 
 nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
                x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
            
[nljac,nljacstr] = symJac(nlcon);
            
 cl = [25;40];
 cu = [Inf;40]; 
 
 %Hessian
 [hess,hstr] = symHessLag(fun,nlcon);

% Bounds (lb <= x <= ub)
 lb = ones(4,1);
 ub = 5*ones(4,1);         

% Initial Guess
 x0 = [1 5 5 1]';
 
opts = optiset('display','iter','solver','ipopt','derivCheck','on','solverOpts',ipoptset('derivative_test','second-order'));
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'H',hess,'HStr',hstr,'bounds',lb,ub,'x0',x0,'options',opts)

[x,f,e,i] = solve(Opt)
 

%% Application to an NLP, analytical derivs supplied
clc
% Objective
 fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Gradient
 grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3)), x(1)*x(4),...
              x(1)*x(4) + 1,  x(1)*sum(x(1:3))];
 
% Nonlinear Constraints 
 nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
                x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
            
% Jacobian
 jac = @(x) sparse([prod(x')./x';
                    2.*x']);

% Jacobian Structure
 jacstr = @() sparse(ones(2,4));
            
 cl = [25;40];
 cu = [Inf;40]; 
 
 %Hessian
% Hessian
 H = @(x,sigma,lambda) sparse(sigma*[ 2*x(4)             0      0   0;
                                      x(4)               0      0   0;
                                      x(4)               0      0   0;
                                      2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ] + ...
                      lambda(1)*[   0          0         0         0;
                                 x(3)*x(4)     0         0         0;
                                 x(2)*x(4) x(1)*x(4)     0         0;
                                 x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ] + ...
                      lambda(2)*diag([2 2 2 2]));

% Hessian Structure
 Hstr = @() sparse(tril(ones(4)));

% Bounds (lb <= x <= ub)
 lb = ones(4,1);
 ub = 5*ones(4,1);         

% Initial Guess
 x0 = [1 5 5 1]';
 
opts = optiset('display','iter','solver','ipopt','derivCheck','on','solverOpts',ipoptset('derivative_test','second-order'));
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'nljac',nljac,'nljacstr',nljacstr,'H',H,'HStr',hstr,'bounds',lb,ub,'x0',x0,'options',opts)

[x,f,e,i] = solve(Opt)
 
 
