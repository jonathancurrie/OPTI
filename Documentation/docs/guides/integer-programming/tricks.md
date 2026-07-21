---
title: "Simple Integer Programming Tricks"
slug: "/guides/integer-programming/tricks/"
---

The following tricks are primarily based on the excellent AIMMS modelling guide:

[https://download.aimms.com/aimms/download/manuals/AIMMS3OM_IntegerProgrammingTricks.pdf](https://download.aimms.com/aimms/download/manuals/AIMMS3OM_IntegerProgrammingTricks.pdf) 

## Constraining a Variable Within a Set
Let's say you have a problem where the variable *x* is an integer (or continuous, it doesn't matter) variable, however *x* may only be certain integer variables within a predefined set. For example:

![int set](/img/opti/int_set.png)

To solve this problem we are going to introduce a binary variable for each valid value within the set. These are known as indicator variables, and when equal to one, indicate the corresponding value in the set is optimal. We will also need a couple of equality constraints to limit the possible values within the set, as well as ensure only one indicator variable is active at the solution.

To illustrate, consider the following simple MILP:

![int set ex](/img/opti/int_set_ex.png)

where we are using the variable *l* as a testing constraint to check the valid values of *x*

```matlab
% Test Constraint Lower Bound
l = 3.5;

% Objective (remember we are not minimizing the indicator vars)
f = [1;0;0;0;0];

% Set Equality Constraints
Aeq = [-1,3,5,8,20    %ensure x in {3,5,8,20}
        0,1,1,1,1];   %only one indicator active
beq = [0;1];

% Bounds
lb = [l;0;0;0;0];
ub = [20;1;1;1;1];

% Integer Declaration
xtype = 'IBBBB';

% Build OPTI Object
Opt = opti('f',f,'eq',Aeq, beq,'bounds',lb,ub,'xtype',xtype)

% Solve
x = solve(Opt)
```

The obvious solution to this problem is the smallest value of *x* that satisfies the bound constraint. However with the addition of the set constraint, it will be the smallest value that is *also* within the set. Try increasing and decreasing *l* to see if the solution always falls within the set.
