---
title: "Mixed Integer Linear Program (MILP)"
slug: "/examples/problem-types/milp/"
---

## Problem Definition
An MILP has the following form:

![def milp](/img/opti/def_milp.png)

Where **f** is a *n* x 1 vector containing the linear objective function, which is subject to the following constraints: 

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

Note an MILP is created in the same way as an LP, except some variables are defined to take discrete (integer) values only. It is recommened you complete reading the [LP](./lp.md) section before reading the remainder of the section.

## Example 1: Small Dense MILP
Consider the following LP from the basics section, modified as an MILP:

![ex1 milp](/img/opti/ex1_milp.png)

Using the native matrix & vector notation of MATLAB this can be entered as so:

```matlab
% Objective
f = -[6 5]';                %Objective Function Vector (min f'x)

% Constraints
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Integer Constraints
xtype = 'II';               %x1 & x2 are Integer
```

And the problem solved by passing the problem variables to OPTI, and calling solve on the resulting object:

```matlab
% Create OPTI Object
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'xtype',xtype)

% Solve the MILP problem
[x,fval,exitflag,info] = solve(Opt)
```

Because the problem contains only two variables, we can use OPTI's built in `plot` command to view the solution:

```matlab
plot(Opt)
```

![plot exmilp](/img/opti/plot_exmilp.png)

You can see from the above figure that the linear gradient of the objective (dashed in grey) as well as the inequality constraints. Shaded areas represent infeasible regions. Unlike an LP, the optimum may not lie on a constraint, where the blue dots indicate feasible integer solutions. Note the Mixed Integer (MI) optimum is different to the relaxed (LP) optimum from the previous example, which is generally the case.

## Example 2: Alternative Integer Setup
There are two ways integer variables can be declared to the OPTI object constructor. 

### Character String
Supply a string of characters representing C (continuous), B (binary) and I (integer):
```matlab
xtype = 'CBI';
```

### Index Vector
Supply a double vector of the indicies of the integer variables:
```matlab
xtype = [2 3];
```

Note that with this method, only integer variables can be defined (not binary). However you can manually add bounds to enforce a binary variable, if required.

## Example 3: Sparse MILP
As with LPs the Example 1 is a toy problem, and real MILP problems will be orders of magnitude bigger. In order to efficiently solve these problems the sparsity structure of the inequality and equality constraints must be exploited.

To illustrate, load a MILP supplied as a .MPS problem with OPTI and lets examine the A matrix:

```matlab
% Load the MILP from .mps file
prob = coinRead('testMILP2.mps')

% Examine A matrix
spy(prob.A)
```

![plot exmilpspy](/img/opti/plot_exmilpspy.png)

While not as sparse as the previous LP example, this problem still contains less than 19% nonzero entries.

To solve this problem, pass the loaded problem structure to the OPTI constructor:

```matlab
% Create an OPTI object from the problem structure
Opt = opti(prob)

% Solve the resulting model
[x,fval] = solve(Opt)
```

## Example 4: MILP with Special Ordered Set (SOS) {#milp-sos}
Consider the following MILP problem:

![ex2 milp](/img/opti/ex2_milp.png)

Initially this problem looks like an LP, however we are going to add one more constraint. We are going to specify that only one variable from *x<sub>1</sub>* to *x<sub>5</sub>* can take on a nonzero value. This type of constraint could be used when choosing one option from a selection.

To model this constraint, we could either introduce binary variables and add the required constraints, or we can specify it as a SOS of type 1. For more information on SOS, consult the [SOS](../../guides/integer-programming/sos.md) section.

Note from OPTI v1.81 the names of SOS structure fields has changed. However you can still enter individual arguments, as is done below:

```matlab
% Objective (f'x)
f = [-1 -1 -3 -2 -2]';

% Row Constraints (rl <= A*x <= ru)
A = sparse([-1 -1  1  1  0;     %sparse A
             1  0  1 -3  0]);
rl = [-Inf;-Inf];    
ru = [30;30];

% Bounds (lb <= x <= ub)
lb = [0;0;0;0;0];
ub = [40;1;Inf;Inf;1];        

% SOS
sos_type = '1'; 
sos_index = [1 2 3 4 5]';
sos_weight = [1 2 3 4 5]';

% Setup Options
opts = optiset('solver','cbc'); %CBC is a SOS constraint solver

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'sos',sos_type,sos_index,...
           sos_weight,'options',opts)

% Solve Problem
[x,fval] = solve(Opt)
```

You will have noticed the introduction of a SOS told OPTI this was a MILP. However we did not need to add any binary variables (or variable types for that matter), and our problem description remains succinct. SOS are also typically handled more efficiently by solvers than if you were to manually include binary variables and constraints, and you can add as many sets as you want!

To add multiple SOS constraints, follow the below construction:
```matlab
% SOS Structure
sos.type = '12'; %each character represents a SOS set 
sos.index = {[1 2]' [3:5]'};
sos.weight = {[1 2]' [1:3]'};

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'sos',sos,'options',opts)
```

Note using the above construct we are utilizing the new SOS API, whereby a single structure is supplied with the fields 'type', 'index', and 'weight'. This API can be used for both single and multiple SOSs, noting cell arrays are used for indices and weights when specifying multiple SOSs. In the above example, we are effectively telling the optimizer to choose one variable between index 1 and 2, and two variables between index 3 and 5.

## Example 5: Specifying Long Integer Variable Strings {#longint}
A common question we get is how to specify `xtype` when you have lots of integer variables. Assuming your variables are ordered (i.e. continuous, integer and binary variables are in consecutive groups), the following example shows a shorthand trick to enter them.

```matlab
% Objective
nC = 10; %Number of Continuous Variables
nI = 10; %Number of Integer Variables
nB = 10; %Number of Binary Variables

% Build xtype vector
xtype = [repmat('C',1,nC),repmat('I',1,nI),repmat('B',1,nB)]
```

## Summary
While the algorithms for solving MILPs have advanced tremendously in the past couple of decades, these are still complex problems and will take much more time than a LP to solve. They also scale poorly with size, so using a commercial solver such as CPLEX which exploits multicore processors can be a significant advantage.

Also note the solution obtained from a MILP solver may not be the global minimum as integer constraints make the problem non-convex.
