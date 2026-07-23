---
title: "Supplying 2nd Derivatives"
slug: "/guides/advanced/deriv2/"
---

Once you are into exact second derivatives you are getting serious! But don't worry, they are not too difficult. However I do suggest reading *all* of the section on [first derivatives](./deriv1.md) before reading this section.

This section will detail how to obtain exact second derivatives using OPTI and MATLAB, and how to supply them to the solver.

## Mathematical Definition
### Scalar Function
The second derivatives of a scalar, multivariable function `f` (such as that used in a NLP objective) is defined as:

$$
\nabla^2 f = \frac{\partial^2 f}{\partial \mathbf{x}^2} = {\def\arraystretch{1.3}\begin{bmatrix} \frac{\partial^2 f}{\partial x_1^2}\; & \;\frac{\partial^2 f}{\partial x_1 \partial x_2}\; & \;\cdots\; & \;\frac{\partial^2 f}{\partial x_1 \partial x_n} \\ \frac{\partial^2 f}{\partial x_2 \partial x_1}\; & \;\frac{\partial^2 f}{\partial x_2^2}\; & \;\cdots\; & \;\frac{\partial^2 f}{\partial x_2 \partial x_n} \\ \vdots\; & \;\vdots\; & \;\ddots\; & \;\vdots \\ \frac{\partial^2 f}{\partial x_n \partial x_1}\; & \;\frac{\partial^2 f}{\partial x_n \partial x_2}\; & \;\cdots\; & \;\frac{\partial^2 f}{\partial x_n^2} \end{bmatrix}}
$$

Noting the result is a *matrix*, with dimensions *n* x *n*.

### Vector Function
The second derivatives of a vector, multivariable function `F` (such as a vector of nonlinear constraints) is a series of Hessian matrices, one for each row, as described above.

### Hessian of the Lagrangian
In practice only one matrix containing second derivative information is passed to a solver. This matrix is called the Hessian of the Lagrangian, and contains second derivatives of *both* the objective *and* constraints:

$$
\nabla^2 L = \sigma \nabla^2 f + \sum_i \lambda_i \nabla^2 c_i
$$

Noting the result has dimensions *n* x *n*.

In the above equation, *f* is the objective function and *c* is the vector of constraints. Two extra arguments are introduced by the solver, sigma, a scalar used to scale the Hessian of the objective, and lambda, a vector where each element scales the Hessian of the respective constraint. The resulting scaled Hessians are summed together to form the Hessian of the Lagrangian.

Note going forward I will refer to the Hessian of the Lagrangian as the Hessian for simplicity.

## Obtaining Exact Second Derivatives
Unfortunately most of the methods listed in the first derivatives section will not work. `autoJac` is not designed for second derivatives, while `mklJac` will be slow and inaccurate. However from OPTI v2.00 you can use `symHess` or `symHessLag`. Both are small 'toy' functions for generating second derivatives (similar to `symJac`), but can be useful for standardised problems. In addition, from OPTI v2.05 you can use `cstepHess` or `cstepHessLag` which implement a modified complex step scheme to estimate second derivatives. These methods are not accurate to numerical precision and can be quite slow, however for small, compatible problems, they work quite well. 

If the methods above do not work alternative options for obtaining 2nd derivatives are:

- Obtain analytical derivatives manually using the Symbolic Toolbox
- A little calculus and enter them manually
- Use [SymBuilder](./sym-builder.md) (under development)
- Consider a 3rd party tool such as  [ADiMat](https://adimat.de/).

However do not fret, IPOPT, BONMIN, and other solvers which use second derivative information are all equipped with a  [Quasi-Newton](http://en.wikipedia.org/wiki/Quasi-Newton_methods) algorithm for approximating the Hessian of the Lagrangian. OPTI enables this functionality automatically if second derivatives are not supplied.

Also note for a particularly dense Hessian, a limited memory  [L-BFGS](http://en.wikipedia.org/wiki/L-BFGS) update which is used to approximate the Hessian can be more efficient that explicitly calculating the elements. This strategy is used by IPOPT and BONMIN. However the inverse is also true, thus for particularly sparse systems it can be quite inefficient.

## Example 1: Hessian Callback Function
Let's revisit Hock & Schittkowski problem #71:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & x_1x_4(x_1+x_2+x_3)+x_3 \\ \text{subject to:} \quad & x_1x_2x_3x_4 \geq 25 \\ & x_1^2+x_2^2+x_3^2+x_4^2=40 \\ & 1 \leq \mathbf{x} \leq 5 \end{aligned}
$$

Assuming we have calculated the Hessian of Lagrangian, it can be entered into an m-file function as follows:

```matlab
function H = hessian (x, sigma, lambda)
%Hessian of the Lagrangian
H = sigma*[ 2*x(4)             0      0   0;
            x(4)               0      0   0;
            x(4)               0      0   0;
            2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ];
H = H + lambda(1)*[    0          0         0         0;
                    x(3)*x(4)     0         0         0;
                    x(2)*x(4) x(1)*x(4)     0         0;
                    x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ];
H = sparse(H + lambda(2)*diag([2 2 2 2]));
```

Note the Hessian has been entered in *lower triangular* form (tril). This the IPOPT and BONMIN default, as it is expected the Hessian is symmetric. MATLAB's `fmincon` on the other hand expects a full Hessian, and OPTI will not convert between the two. Therefore you must unfortunately use the form based on the solver you are using. 

## Example 2: Application to a NLP
Let's supply the entire problem, including all derivatives to OPTI:

```matlab
% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);
 
% Gradient
grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3)), x(1)*x(4),...
             x(1)*x(4) + 1,  x(1)*sum(x(1:3))];

% Nonlinear Constraints 
nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
               x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
cl = [25;40];
cu = [Inf;40]; 

% Jacobian
jac = @(x) sparse([prod(x')./x';
                   2.*x']);

% Jacobian Structure
jacstr = @() sparse(ones(2,4));

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

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
            'hess',H,'hstr',Hstr,'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)
```

Note the Hessian also requires a Hessian structure callback, in exactly the same way the Jacobian structure callback is required. Just remember it should be tril for IPOPT and BONMIN!

## Example 3: Checking 2nd Derivatives
OPTI can also check second derivatives, as follows:

```matlab
% Options
opts = optiset('solver','ipopt','display','iter','derivCheck','on');
 
% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
            'hess',H,'hstr',Hstr,'bounds',lb,ub,'x0',x0,'options',opts)
 
% Solve NLP
[x,fval,exitflag,info] = solve(Opt)
```

If you examine the command window print out, you will see it says (for the above problem):

`OPTI Derivative Checker detected no problems in 'Objective Gradient'`  
`OPTI Derivative Checker detected no problems in 'Constraint Jacobian'`  
`OPTI Derivative Checker detected no problems in 'Hessian of the Lagrangian'`

Indicating no errors were detected.

## Summary
Generating second derivatives by hand is tedious and error prone. However in my experience, especially when solving MINLPs, it can substantially reduce computation time. While the Quasi-Newton approach works well, it is an approximation and often more iterations are required by a solver to find the optimum. Multiply these extra iterations by the number of nodes and relaxed problems to search in a MINLP, and it can quickly add up.
