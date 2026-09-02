---
title: "Quadratically Constrained Quadratic Program (QCQP)"
slug: "/examples/problem-types/qcqp/"
---

## Problem Definition
A QCQP has the following form:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & \tfrac{1}{2}\mathbf{x}^{T}\mathbf{H}\mathbf{x} + \mathbf{f}^{T}\mathbf{x} \\ \text{subject to:} \quad & \mathbf{A}\mathbf{x} \leq \mathbf{b} \\ & \mathbf{A}_{\mathrm{eq}}\mathbf{x} = \mathbf{b}_{\mathrm{eq}} \\ & \mathbf{l}_{\mathrm{b}} \leq \mathbf{x} \leq \mathbf{u}_{\mathrm{b}} \\ & \mathbf{x}^{T}\mathbf{Q}\mathbf{x} + \mathbf{l}^{T}\mathbf{x} \leq r \end{aligned}
$$

Where **H** is a *n* x *n* sparse matrix (quadratic and bilinear terms) and **f** is a *n* x 1 vector (linear terms) containing the quadratic objective function, which is subject to the following constraints: 

### Linear Inequalities*
**A** is a *m* x *n* sparse matrix, **b** is a *m* x 1 vector. 

### Linear Equalities*
**A<sub>eq</sub>** is a *k* x *n* sparse matrix, **b<sub>eq</sub>** is a *k* x 1 vector. 

### Decision Variable Bounds
**l<sub>b</sub>** and **u<sub>b</sub>** are *n* x 1 vectors, where -inf or inf indicate an unbounded lower or upper bound, respectively. 

### Quadratic Constraints
**Q** is a *n* x *n* sparse matrix, **l** is a *n* x 1 vector and *r* is a 1 x 1 scalar. Multiple quadratic constraints are specified by multiple sets of these three variables. **Q** must be convex for all solvers other than SCIP.

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints. 

<small>\*Your problem description will either use Linear Inequalties and Linear Equalities OR Linear Row Constraints. See the [constraint information](../../guides/advanced/cons.md) page.</small>

Note a QCQP is created in a similar way as a QP, so it is recommened you complete reading the [QP](./qp.md) section before reading the remainder of the section.

## Example 1: Small Dense QCQP
Consider the following small QCQP:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & 0.5x_1^2+0.5x_2^2-2x_1-2x_2 \\ \text{subject to:} \quad & -x_1+x_2 \leq 2 \\ & x_1+3x_2 \leq 5 \\ & x_1^2+x_2^2-2x_2 \leq 1 \\ & 0 \leq \mathbf{x} \end{aligned}
$$

Using the native matrix & vector notation of MATLAB this can be entered as so:

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraint
Q = [1 0; 0 1];             %Quadratic Inequality (x'Qx + l'x <= r)
l = [0;-2];
r = 1;
```

And the problem is solved by passing the problem variables to OPTI, and calling solve on the resulting object:

```matlab
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)
```

Because the problem contains only two variables, we can use OPTI's built in `plot` command to view the solution:

```matlab
plot(Opt)
```

![plot exqcqp](/img/opti/plot_exqcqp.png)

The quadratic constraint is the small black circle near the origin, with the hashing indicating the infeasible side. The optimum must always be inside the quadratic constraint circle. As stated, quadratic constraints must be convex (positive definite) unless you are using a global solver, such as [SCIP](../../solvers/scip.md). For more information on convexity, see the [QP page](./qp.md#qpconvex).

## Example 2: QCQP with Multiple Quadratic Constraints {#qcqp-mul}
Supplying multiple quadratic constraints requires the user to stack the individual constraints as follows:

- Q: A cell array of double matrices, column orientated. Each cell is a constraint Q.

- l: A double matrix, where each column is a constraint l vector.

- r: A double column vector, where each row is a constraint r scalar.

To illustrate, consider the following example:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & 0.5x_1^2+0.5x_2^2-2x_1-2x_2 \\ \text{subject to:} \quad & -x_1+x_2 \leq 2 \\ & x_1+3x_2 \leq 5 \\ & x_1^2+x_2^2-2x_2 \leq 1 \\ & x_1^2+x_2^2-x_1+2x_2 \leq 1.2 \\ & 0 \leq \mathbf{x} \end{aligned}
$$

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Inequalities (x'Qx + l'x <= r)
     [1 0; 0 1]};
l = [[0;-2] [-1;2]];
r = [1;1.2];

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)
```

![plot ex2qcqp](/img/opti/plot_ex2qcqp.png)

## Example 3: QCQP with Multiple Quadratic Constraints [Alternate Setup]
Starting from OPTI v1.79 the user can also supply multiple quadratic constraints using the following cell-based format:

- Q: A cell array of double matrices, column orientated. Each cell is a constraint Q.

- l: A cell array of double vectors. Each cell is a constraint l vector.

- r: A cell array of double scalars. Each cell is a constraint r scalar.

To illustrate, we will revisit the above example:

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Inequalities (x'Qx + l'x <= r)
     [1 0; 0 1]};
l = {[0;-2]; [-1;2]};
r = {1; 1.2};

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)
```

## Example 4: Solving a Positive Semidefinite (PSD) QCQP as a NLP {#qcqp-nlp}
OPTI will automatically convert (MI)QCQPs to (MI)NLPs if a compatible NLP solver is requested. Currently [IPOPT](../../solvers/ipopt.md) and [BONMIN](../../solvers/bonmin.md) are setup to solve these problems. The OPTI conversion will automatically generate all first and second derivatives, thus the NLP solver can be quite efficient. 

To illustrate, consider the following example with a PSD quadratic constraint:

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;3.8];    
lb = [0;0];                 %Bounds on x (lb <= x)
ub = [40;inf];

% Quadratic Constraints
Q = [0 0; 0 10];            %Quadratic Inequalities (x'Qx + l'x <= r)
l = [0;-2];
r = 3

% Set OPTI Options
opts = optiset('solver','ipopt','display','iter');

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'bounds',lb,ub,'qc',Q,l,r,'options',opts)

% Solve the QCQP problem [as an NLP]
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt,3)
```

![plot ex4qcqp](/img/opti/plot_ex4qcqp.png)

Note solving Positive SemiDefinite (PSD) QPs or QCQPs does not require you to use a NLP solver. However if your problem is indefinite CPLEX will only solve problems with indefinite terms in the objective (i.e. all quadratic constraints must be PSD). For indefinite quadratic constraints then you must use either [SCIP](../../solvers/scip.md) (a global solver) or use a NLP solver, as above (noting you may only find a local solution using IPOPT or BONMIN).

## Example 5: Quadratic Row Constraints {#qcqp-row}
Starting from OPTI v1.81 you can now supply quadratic constraints in row format, using the following definition:

$$
\mathbf{q}_{\mathrm{rl}} \leq \mathbf{x}^{T}\mathbf{Q}\mathbf{x} + \mathbf{l}^{T}\mathbf{x} \leq \mathbf{q}_{\mathrm{ru}}
$$

Note the above format adds extra reasons why the constraint may not be convex. This could be due to:

1. If qrl = qru, then the quadratic constraint is an equality, and nonlinear equality constraints are non-convex.
1. If both qrl and qru are finite (i.e. a double-sided constraint), then by definition one side will result in a concave inequality.
1. As per standard quadratic inequalities, if Q is not positive semi-definite, then the constraint is not convex.

If any of the above conditions are true then you are now solving a non-convex QCQP, which is considerably more difficult to solve! However if you have [SCIP](../../solvers/scip.md) installed (a global solver), OPTI will automatically use it to solve your model, as shown below:

```matlab
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Constraints (qrl <= x'Qx + l'x <= qru)
     [1 0; 0 1]};
l = {[0;-2]; [-2;2]};
qrl = {3; 1};               %QC1 is double sided, QC2 is an equality
qru = {5; 1};

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qcrow',Q,l,qrl,qru)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)
```

![plot ex5qcqp](/img/opti/plot_ex5qcqp.png)

As seen in the above plot, the double-sided quadratic constraint forms a 'donut' shaped feasible region. The quadratic equality is shown in blue, as per the default nonlinear equality constraint colours.

## OPTI QCQP Solvers
If you are an academic user and have downloaded the academic version of OPTI a new solver, [SCIP](../../solvers/scip.md), is now included which can solve convex and non-convex QCQPs. Alternatively you can download CPLEX which can solve convex QCQPs to global optimality.

If you do not have these installed OPTI will automatically convert the problem to a MINLP and solve it with IPOPT if no dedicated QCQP solver is available. While not as efficient, all first and second derivative information is passed to IPOPT and it does quite well at solving these problems.

## Summary
Solving QCQPs is a complex problem thus long computation times can be expected on large scale problems. Note the 0.5 is dropped from the quadratic constraint definition in OPTI (CPLEX convention).
