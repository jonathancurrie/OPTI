%% OPTI Toolbox Differentiation Examples
%
% This file demonstrates each of the different differentiation algorithms
% supplied with the OPTI Toolbox. Note you should be familiar with the 
% operation of OPTI Toolbox by reading the accompanying examples.
%
%   Copyright (C) 2014 Jonathan Currie (Control Engineering)

%Also see the following webpage for more examples:
web('https://www.controlengineering.co.nz/Wikis/OPTI/index.php/Advanced/Deriv1');

%% Gradient vs Jacobian
% The terms Gradient and Jacobian are associated with differentiation
% functions, and for the purposes of this toolbox they are basically
% interchangeable:
%
% Gradient - The first derivative of the objective function (1 x n)
% Jacobian - The first derivative of the constraint function (m x n)
%
% Where you will note the main difference is the gradient is a vector, and
% the Jacobian is a matrix. However this is just a rule of thumb and will
% not hold in all instances (constraint functions with one row for
% example!). All functions are named as returning a Jacobian, which could
% also be a Gradient Vector.
 
%% Problem 1
% For the following examples we are going to use this nonlinear objective
% function:
clc
fun = @(x) sin(pi*x(1)/12) * cos(pi*x(2)/16);

x = [0.75 0.25]';

%% Example 1 - Automatic Differentiation (AD)
% Automatic differentiation is one of the most powerful differentiation
% strategies which can provide error free gradients of a function by
% applying the chain rule to each operation. 
%
% The AD routines implemented in OPTI are supplied by adiff, a Matlab
% project by William McIlhagga. While most Matlab functions are overloaded,
% you cannot use external code (i.e. via MEX) or any toolbox or class
% functions.

dx = autoJac(fun,x)

%% Example 2 - Numerical Differentiation (ND) [Finite Difference]
% Numerical differentiation using finite differences is a computationally
% expensive procedure which can result in an inaccurate gradient if the
% internal perturbations are not chosen correctly. 
%
% The ND routine implemented in OPTI is the Intel MKL djacobi function
% which approximates the derivative using central differences. This is
% implemented via a MEX function which repeatedly calls the function in
% order to close in on the gradient.

dx = mklJac(fun,x)

%% Example 3 - Numerical Differentiation (ND) [Complex Step]
% Complex step differentiation works on native MATLAB functions, provided
% you have your ' and .' not mixed up (conjugate vs transpose!)

dx = cstepJac(fun,x)

%% Example 4 - Symbolic Differentiation (SD)
% Symbolic differentiation analytically differentiates the function as a
% symbolic expression, resulting in a single expression for the gradient.
% Complications occur if the function cannot be analytically differentiated
% or the symbolic routine cannot find a derivative.
%
% The SD routine implemented in OPTI uses the Matlab Symbolic Toolbox as
% well as two wrapper functions in order to generate the gradient. The
% wrapper functions convert the function handle to a symbolic expression
% and vice-versa.

grad = symJac(fun)

if(~isempty(grad)) %don't run if Symbolic Toolbox not installed
    dx = grad(x)
end


%% Example 5 - Applying AD to NLP Solving
% To use AD to generate the objective gradient for IPOPT create a function
% handle which calls autoJac, and then pass this to the optiprob function:
clc
%Problem
obj = @(x) log(1+x(1)^2) - x(2);
lb = [-2 -2]';
ub = [2 2]';

grad = @(x) autoJac(obj,x);

opts = optiset('solver','ipopt');
Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

[x,fval,exitflag,info] = solve(Opt,[0;0]);
fval
info

%% Example 6 - Applying ND to NLP Solving
% To use ND to generate the objective gradient for IPOPT just replace the
% above with mklJac:
clc
grad = @(x) mklJac(obj,x);

Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

[x,fval,exitflag,info] = solve(Opt,[0;0]);
fval
info

%% Example 6 - Applying ND to NLP Solving
% To use ND to generate the objective gradient for IPOPT just replace the
% above with cstepJac:
clc
grad = @(x) cstepJac(obj,x);

Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

[x,fval,exitflag,info] = solve(Opt,[0;0]);
fval
info

%% Example 7 - Applying SD to NLP Solving
% To use SD we can generate the analytical gradient and use this as our
% gradient function:
clc
grad = symJac(obj)

if(~isempty(grad))
    Opt = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',opts)

    [x,fval,exitflag,info] = solve(Opt,[0;0]);
    fval
    info
else
    display('Cannot run example without Symbolic Toolbox!');
end

%% Example 8 - Ensuring Accurate Derivatives
% OPTI comes with a built in derivative checker so you can check your
% derivatives. Enable it via the 'derivCheck' option in optiset. NOTE
% supply x0 to opti() to check the derivatives at the user supplied start
% point.

% Includes Hessian
clc
%Problem
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)]';
          
%Hessian of the Lagrangian & Structure       
hess = @(x,sigma,lambda) sparse(sigma*[ 1200*x(1)^2-400*x(2)+2  0          0                          0
                                        -400*x(1)               220.2      0                          0
                                         0                      0          1080*x(3)^2-360*x(4)+2     0
                                         0                      19.8       -360*x(3)                  200.2 ]);  
Hstr = @() sparse([ 1  0  0  0 
                    1  1  0  0
                    0  0  1  0
                    0  1  1  1 ]);   
                
%Bounds
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
x0 = [-3  -1  -3  -1]';

% The Hessian + Structure are added just like other arguments. However note
% some solvers require a TRIL, TRIU or SYM Hessian. For instance IPOPT uses
% TRIL (assumes symmetric), while FMINCON requires the full Hessian. Make
% sure to configure your inputs as per the solver you're using.
opts = optiset('derivCheck','on');
Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub,'x0',x0,'options',opts)

[x,fval,exitflag,info] = solve(Opt)

%% Conclusion
% By default OPTI Toolbox will always use Numerical Differentiation
% (mklJac) for gradients / Jacobians if one has not been provided. This
% provides the best combination of flexibility with complex objective / 
% constraint functions and performance. However if you can use SD to
% generate a gradient this would be a preferred method. 
%
% If you find the optimizer is failing with ND and AD is a suitable 
% candidate you can try it to ensure your gradients are accurate.
