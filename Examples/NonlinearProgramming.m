%% OPTI Toolbox Nonlinear Programming Examples
%
% This file contains a number of NLP problems and demonstrates how to 
% solve them using the OPTI Toolbox. You should read and complete
% BasicUsage.m & LinearProgramming.m BEFORE running the below examples.
%
%   Copyright (C) 2014 Jonathan Currie (IPL)

% There is also a page on the Wiki which supplements this example:
web('https://www.inverseproblem.co.nz/OPTI/index.php/Probs/NLP');

%% Determing which Solver to Use
% OPTI Toolbox comes with a number of NLP solvers, thus to determine which
% ones are available on your system you can type:
clc
optiSolver('NLP')

%% NLP Considerations - 1st Derivatives
% See DifferentationExamples.m

%% NLP Considerations - Helping the Solver
% The most important point when solving NLPs is to give the solver as much
% information as you know, this should not only speed up the solver, but
% may also result in a better optimum. OPTI Toolbox allows the following
% (but not all solvers may support all the below):

% - Gradient Function (f,grad)
% A function handle which the optimizer can use to determine the 1st
% derivatives of the objective function.

% - Hessian Function (H,hess)
% A function handle which the optimizer can use to determine the 2nd
% derivatives of the objective function + constraint functions.

% - Hessian Structure (Hstr)
% A function handle of returning a sparse matrix indicating the non-zero
% terms within the Hessian. This saves a numerical difference / BFGS update
% from evaluating these terms in the Hessian.

% - Jacobian Function (nljac)
% A function handle which the optimizer can use to determine the 1st
% derivatives of the constraint function.

% - Jacobian Structure (nljacstr)
% A function handle of returning a sparse matrix indicating the non-zero
% terms within the Jacobian. 

% And most importantly - use a good (i.e. feasible & well chosen) initial
% guess for x0!

%% Example 1
% This is a simple two decision variable NLP which will use for the next 
% few examples
clc

%Problem
obj = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2; %nl objective
grad = @(x)[2*x(1)-400*x(1)*(x(2)-x(1)^2)-2;     %nl gradient
            200*x(2)-200*x(1)^2];

lb = [-inf; -1.5];      %Bounds on x (lb <= x <= ub) 
x0 = [-2; 1];           %Starting guess of solution               

% Build OPTI Problem
Opt = opti('obj',obj,'grad',grad,'lb',lb)

% Call solve to solve the problem
[x,fval,exitflag,info] = solve(Opt,x0)   

%% Example 2 - Alternative Setup Strategies
% Naming of arguments, as well as pairing is flexible when using opti

Opt = opti('fun',obj,'grad',grad,'lb',lb); %fun = obj

% OR
Opt = opti('obj',obj,'f',grad,'lb',lb); %f = grad

%% Example 3 - Using Numerical Difference for Gradient
% Simply drop 'grad' from the problem to use a numerical difference scheme
% instead - you will note mklJac (a numerical difference routine) is now 
% used for the gradient.
clc
Opt = opti('obj',obj,'lb',lb)

x = solve(Opt,x0)  

%% Example 4 - Unconstrained Nonlinear Optimization
% If you do not supply any constraints OPTI will solve the problem using an
% unconstrained NL solver. For problems where OPTI cannot determine the
% number of decision variables you must supply this to optiprob:
clc
Opt = opti('obj',obj,'ndec',2)

x = solve(Opt,x0)  

%% Example 5 - Nonlinear Constraints
% Includes Nonlinear Constraints
clc
%Problem
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
            
%Nonlinear Constraints + Jacobian + Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) sparse([ 1  3  0  0  0;
	                  0  0  1  1 -2;
	                  0  1  0  0 -1 ]);
nljacstr = @() sparse([1 1 0 0 0;
                       0 0 1 1 1;
                       0 1 0 0 1]);
nlrhs = [4 0 0]';   %nonlinear constraints RHS
nle = [0 0 0]';     %nonlinear constraints type (-1 <=, 0 ==, 1 >=)
x0 = [ 2.5 0.5 2 -1 0.5 ];

% Nonlinear constraints are added using 'nlmix' or individual args
% ('nlcon','nlrhs','nle'). To remain flexible the right hand side of the
% equations is set via nlrhs (some solvers always default to 0) as well as
% the constraint type (>=,<=,==) via nle.
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr)

[x,fval,exitflag,info] = solve(Opt,x0)  

%% Example 6 - User Supplied Hessian
% Includes Hessian
clc
%Problem
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];
          
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
Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub)

[x,fval,exitflag,info] = solve(Opt,x0)

%% A Note on Hessians (2nd Derivatives)
% Be aware the Hessian required by OPTI is the Hessian of the Lagrangian.
% This means it contains the second derivatives of both the objective
% (scaled by sigma) and the nonlinear constraints (scaled by lambda). 
%
% An example Hessian callback including nonlinear constraints (HS #71) 
% is shown below:
%
%      function H = hessian (x, sigma, lambda)
%           H = sigma*[ 2*x(4)             0      0   0;
%                       x(4)               0      0   0;
%                       x(4)               0      0   0;
%                       2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
%           H = H + lambda(1)*[    0          0         0         0;
%                               x(3)*x(4)     0         0         0;
%                               x(2)*x(4) x(1)*x(4)     0         0;
%                               x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
%           H = sparse(H + lambda(2)*diag([2 2 2 2]));

%% Example 7 - No Derivatives
% With no gradient or Jacobian supplied, a numerical difference algorithm
% is used for both. Note ndec must be supplied
clc
%Problem
obj = @(x) sin(pi*x(1)/12)*cos(pi*x(2)/16);
nlcon = @(x) 4*x(1)-3*x(2);
nlrhs = 0;
nle = 0;
x0 = [0;0];

Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'ndec',2)

[x,fval,exitflag,info] = solve(Opt,x0)

%% Example 8 - Derivative Free Optimization
% NLOPT comes with several derivative free optimizers, but not many can
% solve with nonlinear constraints. Use the command below to see all
% algorithms and statistics:
clc
nloptSolver

%% Example 8 - Building an NLOPT Problem
% Lets choose LN_COBYLA to solve the problem and use nloptset to select it:
clc
nopt = nloptset('algorithm','LN_COBYLA');            %nloptset functions much like optiset
opts = optiset('solver','nlopt','solverOpts',nopt);  %specify solver + our NLOPT options above

Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'ndec',2,'options',opts);

[x,fval,exitflag,info] = solve(Opt,x0)

%% Example 9 - Advanced NLOPT Setup
% You can specify advanced NLOPT options such as the Augmented Langragian
% Solver (for constrained problems but solved with a bounded suboptimzer)
% as follows:
clc
%Problem
obj = @(x) x(1) - x(2);
nlcon = @(x) -3*x(1)^2 + 2*x(1)*x(2) - x(2)^2 + 1;
nlrhs = 0;
nle = 1;
x0 = [-10;10];

nopt = nloptset('algorithm','LN_AUGLAG','subalgorithm','LN_NELDERMEAD'); %add the suboptimzer settings        
opts = optiset('solver','nlopt','solverOpts',nopt);  
%Build OPTI Problem
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'x0',x0,'options',opts)
%Solve
[x,fval,exitflag,info] = solve(Opt)

