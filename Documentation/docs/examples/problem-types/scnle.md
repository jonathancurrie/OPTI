---
title: "System of Constrained Nonlinear Equations (SCNLE)"
slug: "/examples/problem-types/scnle/"
---

## Problem Definition
A SCNLE has the following form:

$$
\begin{aligned} & \mathbf{F}(\mathbf{x}) = \mathbf{0} \\ \text{subject to:} \quad & \mathbf{A}\mathbf{x} \leq \mathbf{b} \\ & \mathbf{A}_{\mathrm{eq}}\mathbf{x} = \mathbf{b}_{\mathrm{eq}} \\ & \mathbf{l}_{\mathrm{b}} \leq \mathbf{x} \leq \mathbf{u}_{\mathrm{b}} \end{aligned}
$$

Where **F** is a vector function containing the nonlinear equations.

### Linear Inequalities
**A** is a *m* x *n* dense matrix, **b** is a *m* x 1 vector. 

### Linear Equalities
**A<sub>eq</sub>** is a *k* x *n* dense matrix, **b<sub>eq</sub>** is a *k* x 1 vector. 

### Decision Variable Bounds
**l<sub>b</sub>** and **u<sub>b</sub>** are *n* x 1 vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively. 

The goal is to set the function values of all equations to zero by selecting a value of **x** that also satisfies all constraints. This problem is known as constrained root solving, or multivariable root solving when the dimension of **x** is greater than 1.

Note a SCNLE is created in a similar way as a SNLE problem. It is recommened you complete reading the [SNLE](./snle.md) section before reading the remainder of the section.

## Example 1: Bounded SNLE

$$
\begin{aligned} 2x_1-x_2-e^{-x_1} &= 0 \\ -x_1+2x_2-e^{-x_2} &= 0 \\ \text{subject to:} \quad 0.6 \leq x_1 &\leq 1 \\ 0 \leq x_2 &\leq 1 \end{aligned}
$$

To setup this problem, it can be entered as follows:

```matlab
% System of Nonlinear Equations
nleq = @(x) [ 2*x(1) - x(2) - exp(-x(1));
            -x(1) + 2*x(2) - exp(-x(2))];

% Bounds
lb = [0.6;0];
ub = [1;1];
 
% Starting Guess
x0 = [-5;5];
```

And the problem solved by passing the problem variables to OPTI, and calling solve on the resulting object:

```matlab
% Create OPTI Object
Opt = opti('nleq',nleq,'bounds',lb,ub,'x0',x0)

% Solve the SCNLE problem
[x,fval,exitflag,info] = solve(Opt)
```
