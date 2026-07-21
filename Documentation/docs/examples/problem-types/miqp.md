---
title: "Mixed Integer Quadratic Program (MIQP)"
slug: "/examples/problem-types/miqp/"
---

## Problem Definition
An MIQP has the following form:

![def miqp](/img/opti/def_miqp.png)

Where **H** is a *n* x *n* sparse matrix (quadratic and bilinear terms) and **f** is a *n* x 1 vector (linear terms) containing the quadratic objective function, which is subject to the following constraints: 

### Linear Inequalities*
**A** is a *m* x *n* sparse matrix, **b** is a *m* x 1 vector. 

### Linear Equalities*
**A<sub>eq</sub>** is a *k* x *n* sparse matrix, **b<sub>eq</sub>** is a *k* x 1 vector. 

### Decision Variable Bounds
**l<sub>b</sub>** and **u<sub>b</sub>** are *n* x 1 vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively. 

### Integer Constraints
*x<sub>i</sub>* are decision variables which must be an integer number (..-2, -1, 0, 1, 2..).

### Binary Constraints
*x<sub>j</sub>* are decision variables which must be a binary number (0,1), where i ≠ j.

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints. 

<small>\*Your problem description will either use Linear Inequalties and Linear Equalities OR Linear Row Constraints. See the [constraint information](../../guides/advanced/cons.md) page.</small>

Note a MIQP is created in a very similar way as a QP, so it is recommened you complete reading the [QP](./qp.md) section before reading the remainder of the section.

## Example 1: Small MIQP
Consider the following small MIQP modified from the problem in the QP section:

![ex1 miqp](/img/opti/ex1_miqp.png)

Using the native matrix & vector notation of MATLAB this can be entered as so:

```matlab
% Objective
H = [1 -1; -1  2];          %Objective Function (min 0.5x'Hx + f'x)
f = -[2 6]';                

% Constraints
A = [1,1; -1,2; 2, 1];      %Linear Inequality Constraints (Ax <= b)
b = [2;2;3];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Integer Constraints
xtype = 'IC';
```

And the problem is solved by passing the problem variables to OPTI, and calling solve on the resulting object:

```matlab
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'xtype',xtype)

% Solve the MIQP problem
[x,fval,exitflag,info] = solve(Opt)
```

Because the problem contains only two variables, we can use OPTI's built in `plot` command to view the solution:

```matlab
plot(Opt,2)
```

![plot exmiqp](/img/opti/plot_exmiqp.png)

## OPTI MIQP Solvers
If you are an academic user and have downloaded the academic version of OPTI a new solver, [SCIP](../../solvers/scip.md), is now included which can solve convex and non-convex MIQPs. Alternatively you can download CPLEX which can solve convex MIQPs and non-convex MIQPs to first-order optimality.

If you do not have these installed OPTI will automatically convert the problem to a MINLP and solve it with BONMIN if no dedicated MIQP solver is available. While not as efficient, all first and second derivative information is passed to BONMIN and it does quite well at solving these problems.

## Summary
Just like with MILPs and LPs, solving a MIQP is much more computationally expensive than solving a QP.
