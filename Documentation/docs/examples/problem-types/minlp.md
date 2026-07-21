---
title: "Mixed Integer Nonlinear Program (MINLP)"
slug: "/examples/problem-types/minlp/"
---

## Problem Definition
An MINLP has the following form:

![def minlp](/img/opti/def_minlp.png)

Where *f* is a scalar function containing the nonlinear objective function, which is subject to the following constraints: 

### Linear Inequalities*
**A** is a *m* x *n* sparse matrix, **b** is a *m* x 1 vector. 

### Linear Equalities*
**A<sub>eq</sub>** is a *k* x *n* sparse matrix, **b<sub>eq</sub>** is a *k* x 1 vector. 

### Decision Variable Bounds
**l<sub>b</sub>** and **u<sub>b</sub>** are *n* x 1 vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively. 

### Nonlinear Inequalities*
**c** is a *u* x 1 vector of functions containing nonlinear inequality constraints, **d** is a *u* x 1 vector 

### Nonlinear Equalities*
**c<sub>eq</sub>** is a *v* x 1 vector of functions containing nonlinear equality constraints, **d<sub>eq</sub>** is a *v* x 1 vector

### Integer Constraints
*x<sub>i</sub>* are decision variables which must be an integer number (..-2, -1, 0, 1, 2..).

### Binary Constraints
*x<sub>j</sub>* are decision variables which must be a binary number (0,1), where i ≠ j.

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints. 

<small>\*Your problem description will either use Linear / Nonlinear Inequalties and Linear / Nonlinear Equalities OR Linear / Nonlinear Row Constraints. See the [constraint information](../../guides/advanced/cons.md) page.</small>

Note an MINLP is created in the same way as an NLP, except some variables are defined to take discrete (integer) values only. It is recommened you complete reading the [NLP](./nlp.md) section before reading the remainder of the section.

## Example 1: Convex MINLP
Before you read any further, have another think about the term "Convex MINLP". By introducing integer constraints, the problem is now non-convex. So how can an MINLP be convex? Well, technically it's not. However the relaxed problem (integer constraints removed) is convex, thus this is typically what is meant by "Convex MINLP". Just remember it really is technically incorrect!

Anyway, consider the following MINLP:

![ex1 minlp](/img/opti/ex1_minlp.png)

This problem contains a linear objective, one quadratic constraint and two linear constraints. While it could be formulated as a MIQCQP, we will keep it as a MINLP for this example.

The problem is entered into MATLAB as follows:
```matlab
% Objective
fun = @(x) -x(1) - x(2) - x(3);    

% Linear Constraints
A = [1 -1 0 0;
     1  0 1 1];
b = [0;2];

% Nonlinear Constraint
nlcon = @(x) (x(2) - 0.5)^2 + (x(3) - 0.5)^2; 
nlrhs = 0.25;
nle = -1; % -1 for <=, 0 for ==, +1 >=         

% Bounds
lb = [0;0;0;0];
ub = [1;Inf;Inf;5];

% Integer Constraints
xtype = 'BCCI';

%Initial Guess
x0 = [0;0;0;0];                      
```

Note we have opted to split the constraints into linear and nonlinear types. While the MINLP solver will call them all in a single nonlinear callback, OPTI will recognise two as being linear, and pass this linear information to the solver to use while solving.

To solve the problem:

```matlab
% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'bounds',lb,ub,...
           'xtype',xtype)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)
```

Using BONMIN (the default MINLP solver), this problem should be solved at the root node, requiring less than 1/4 of a second to solve. However remember this is a tiny problem, real MINLPs can take hours or days to solve.

## Example 2: Non-Convex MINLP
Like IPOPT, BONMIN will only find local solutions to non-convex problems, thus for non-convex MINLPs we have to use a different solver. Rastrigin's function is an example of a test non-convex NLP:

![plot ex2minlp](/img/opti/plot_ex2minlp.png)

To make it more interesting, we are going to constrain x<sub>1</sub> to be an integer variable, within a curved region of the function:

![ex2 minlp](/img/opti/ex2_minlp.png)

This can be entered into MATLAB as follows:

```matlab
% Objective
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)))

% Bounds
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];

% Integer Constraints
xtype = 'IC';

%Initial Guess
x0 = [16;0];                      
```

For this problem we are going to use NOMAD to solve it. NOMAD is an excellent derivative free, global MINLP solver which continues to surprise us with it's ability to find solutions to difficult, real world problems.

NOMAD is used to solve the problem as follows:
```matlab
% Options
opts = optiset('solver','nomad','display','iter')

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'xtype',xtype,'options',opts)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)                    
```

NOMAD should return in less than 1/4 of a second with a integer optimal solution. Note we have also enabled iteration display, allowing us to view solver progress while it is solving. On large problems this is especially useful, to ensure progress is being made!
