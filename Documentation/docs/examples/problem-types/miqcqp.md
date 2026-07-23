---
title: "Mixed Integer Quadratically Constrained Quadratic Program (MIQCQP)"
slug: "/examples/problem-types/miqcqp/"
---

## Problem Definition
An MIQCQP has the following form:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & \tfrac{1}{2}\mathbf{x}^{T}\mathbf{H}\mathbf{x} + \mathbf{f}^{T}\mathbf{x} \\ \text{subject to:} \quad & \mathbf{A}\mathbf{x} \leq \mathbf{b} \\ & \mathbf{A}_{\mathrm{eq}}\mathbf{x} = \mathbf{b}_{\mathrm{eq}} \\ & \mathbf{l}_{\mathrm{b}} \leq \mathbf{x} \leq \mathbf{u}_{\mathrm{b}} \\ & \mathbf{x}^{T}\mathbf{Q}\mathbf{x} + \mathbf{l}^{T}\mathbf{x} \leq r \\ & x_i \in \mathbb{Z} \\ & x_j \in \{0,1\} \end{aligned}
$$

Where **H** is a *n* x *n* sparse matrix (quadratic and bilinear terms) and **f** is a *n* x 1 vector (linear terms) containing the quadratic objective function, which is subject to the following constraints: 

### Linear Inequalities*
**A** is a *m* x *n* sparse matrix, **b** is a *m* x 1 vector. 

### Linear Equalities*
**A<sub>eq</sub>** is a *k* x *n* sparse matrix, **b<sub>eq</sub>** is a *k* x 1 vector. 

### Decision Variable Bounds
**l<sub>b</sub>** and **u<sub>b</sub>** are *n* x 1 vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively. 

### Quadratic Constraints
**Q** is a *n* x *n* sparse matrix, **l** is a *n* x 1 vector and **r** is a 1 x 1 scalar. Multiple quadratic constraints are specified by multiple sets of these three variables. **Q** must be convex.

### Integer Constraints
*x<sub>i</sub>* are decision variables which must be an integer number (..-2, -1, 0, 1, 2..).

### Binary Constraints
*x<sub>j</sub>* are decision variables which must be a binary number (0,1), where i ≠ j.

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints. 

<small>\*Your problem description will either use Linear Inequalties and Linear Equalities OR Linear Row Constraints. See the [constraint information](../../guides/advanced/cons.md) page.</small>

Note a MIQCQP is created in a similar way as a QP and QCQP, so it is recommened you complete reading the [QP](./qp.md) and [QCQP](./qcqp.md) sections before reading the remainder of the section.

## Example 1: Small Dense MIQCQP
Consider the following small MIQCQP:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & 0.5x_1^2+0.5x_2^2-2x_1-2x_2 \\ \text{subject to:} \quad & -x_1+x_2 \leq 2 \\ & x_1+3x_2 \leq 5 \\ & x_1^2+x_2^2-2x_2 \leq 1 \\ & 0 \leq \mathbf{x} \\ & x_1 \in \mathbb{Z} \end{aligned}
$$

Using the native matrix & vector notation of MATLAB this can be entered as so:

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraint
Q = [1 0; 0 1];             %Quadratic Inequality (x'Qx + l'x <= r)
l = [0;-2];
r = 1;

% Integer Constraints
xtype = 'IC';
```

And the problem is solved by passing the problem variables to OPTI, and calling solve on the resulting object:

```matlab
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r,'xtype',xtype)

% Solve the MIQCQP problem
[x,fval,exitflag,info] = solve(Opt)
```

Because the problem contains only two variables, we can use OPTI's built in `plot` command to view the solution:

```matlab
plot(Opt,3)
```

![plot exmiqcqp](/img/opti/plot_exmiqcqp.png)

The quadratic constraint is the small black circle near the origin, with the hashing indicating the infeasible side. As only x<sub>1</sub> is required to be an integer, x<sub>2</sub> takes a value outside the blue integer dot.

## OPTI MIQCQP Solvers
If you are an academic user and have downloaded the academic version of OPTI a new solver, [SCIP](../../solvers/scip.md), is now included which can solve convex and non-convex MIQCQPs. Alternatively you can download CPLEX which can solve convex MIQCQPs.

If you do not have these installed OPTI will automatically convert the problem to a MINLP and solve it with BONMIN if no dedicated MIQCQP solver is available. While not as efficient, all first and second derivative information is passed to BONMIN and it does quite well at solving these problems.
