---
title: "MATLAB Optimization Toolbox Overloads"
slug: "/getting-started/overloads/"
---

It is expected many users will be experienced with the MATLAB Optimization Toolbox. Therefore to minimize the transition of code between the two, OPTI provides several *overloads*. Note these are not overloads in the object-orientated sense, merely the same function names with 'opti_' in front of it.

The following table lists the Optimization Toolbox and corresponding OPTI functions:

|  |  |  |
| --- | --- | --- |
| Problem Type | MATLAB Function | OPTI Function |
| LP | `linprog` | `opti_linprog` |
| BILP | `bintprog` | `opti_bintprog` |
| MILP | `intlinprog` | `opti_intlinprog` |
| QP | `quadprog` | `opti_quadprog` |
| NLS | `lsqcurvefit` | `opti_lsqcurvefit` |
| SNLE | `fsolve` | `opti_fsolve` |
| UNO | `fminunc` | `opti_fminunc` |
| NLP | `fmincon` | `opti_fmincon` |

Each function uses exactly the same calling syntax as the MATLAB equivalent, including dual outputs for nonlinear callbacks (not an OPTI convention).

## Example 1: LP
Using the same LP from the Basic Usage example page, the problem can be constructed and solved as follows:

```matlab
% Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Solve using MATLAB's linprog:
x = linprog(f,A,b,[],[],lb,ub)

% Solve using corresponding OPTI overload:
x = opti_linprog(f,A,b,[],[],lb,ub)
```

Note OPTI is not intended to replace the Optimization Toolbox. It in fact interfaces with it, and you can choose 'matlab' as the solver for any problem type you wish to attempt to solve using the Optimization Toolbox. In my experience the Optimization Toolbox is a suite of robust solvers and can typically solve problems others can't. It is just not always that fast, which can be a problem on large-scale problems.

## Example 2: MILP
From MATLAB 2014a you can now solve mixed-integer linear programs. However if you want to try one of OPTI's solvers, use the overload below:
```matlab
% Objective
f = [2, 3, 7, 7]';         
% Constraints       
A = [-1, -1, 2, 5;1, -2, -1, -4];
b = [-2; -3];  
lb = zeros(4,1); 
ub = [30 100 20 1]';  

% Integer Constraints
xtype = [2 4];  %variables 1, 3 are continuous, 2, 4 are integer

% Solve using OPTI overload:
x = opti_intlinprog(f,xtype,A,b,[],[],lb,ub)
```

Note an earlier overload exists, `opti_mintprog`, which contains the same functionality. This function existed before `intlinprog` was released.

## Example 3: NLP
You may be wondering how to change the solver using the overloaded methods. I can assure you this is easy! Just as you can pass optimset options to a MATLAB function, you can pass optiset options (including the solver) to a OPTI overload. The following example details this:

```matlab
% Objective
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;          
% Constraints       
lb = [-1.5;-3];
ub = [4;3];  

% Starting Guess
x0 = [0;0];

% Setup Options
opts = optiset('solver','nomad');

% Solve using OPTI overload:
[x,fval] = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)
```

## Summary
These overloads are provided for convenience only. They do not represent the full functionality of either the Optimization Toolbox or OPTI, just a small subset so you can interchange between them easily for testing purposes. It is strongly recommended you follow standard OPTI syntax for normal use.
