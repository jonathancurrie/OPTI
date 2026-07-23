---
title: "Global Nonlinear Program (GNLP)"
slug: "/examples/problem-types/gnlp/"
---

## Problem Definition
An GNLP has the following form:

$$
\begin{aligned} \min_{\mathbf{x}} \quad & f(\mathbf{x}) \\ \text{subject to:} \quad & \mathbf{A}\mathbf{x} \leq \mathbf{b} \\ & \mathbf{A}_{\mathrm{eq}}\mathbf{x} = \mathbf{b}_{\mathrm{eq}} \\ & \mathbf{l}_{\mathrm{b}} \leq \mathbf{x} \leq \mathbf{u}_{\mathrm{b}} \\ & \mathbf{c}(\mathbf{x}) \leq \mathbf{d} \\ & \mathbf{c}_{\mathrm{eq}}(\mathbf{x}) = \mathbf{d}_{\mathrm{eq}} \end{aligned}
$$

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

The goal is to minimize the objective function by selecting a value of **x** that also satisfies all constraints. 

<small>\*Your problem description will either use Linear / Nonlinear Inequalties and Linear / Nonlinear Equalities OR Linear / Nonlinear Row Constraints. See the [constraint information](../../guides/advanced/cons.md) page.</small>

Note an GNLP is created in the exact same way as an NLP. It is recommened you complete reading the [NLP](./nlp.md) section before reading the remainder of the section.

## GNLP vs NLP
You may have noticed the defintion of a GNLP is exactly the same as a NLP. Well in fact it is (in my terminology), what varies is the how the solver works, and the type of solution returned.

A global solver will attempt to find the *global* optimum to a problem, which may be significantly different from a local solution. There are multiple strategies to do this, however they fall into two general categories: derivative based and derivative free.

### Derivative Based Global Optimizers
Global optimizers which use derivatives must use extra information, other than just the derivatives, in order to find a solution, otherwise they would fall into the same local minima. Typically a white-box model is required, where the model, including all the equations, is passed directly to the solver. This allows the optimizer to exploit the structure of the system in order to find a global solution, typically using a cutting plane based approach together with branch and bound.

As of OPTI v1.73 SCIP has been updated to allow it to solve NLPs and MINLPs to global optimality. The interface parses your MATLAB functions into an algebraic description, which is then used by the solver, see the [White Box Optimization](../../guides/advanced/white-box.md) page for more details. In addition, from OPTI v2.05 the global optimization solver BARON is also interfaced. There is also an open-source project, Couenne, which implements this functionality, but I have not had the time to write an interface for it.

An alternative to the white-box approach is multi-start. This runs a local optimizer over many initial start points, and returns the best solution from all runs. The [multisolve](../../guides/advanced/multi-solve.md) algorithm supplied with OPTI provides a very basic implementation of a multi-start solver.

### Derivative Free Global Optimizers
Derivative free optimizers use strategies like direct search, particle swarm, ant colony, and genetic algorithm to solve the problem. These solvers require no derivative information, and instead only require objective and constraint evaluations. The downside is they will typically use *many* more function evaluations in order to find a solution.

Derivative free optimizers are also suited to general NLP problems, where the objective or constraints may be non-differentiable, or subject to stochastic noise. An example is fitting parameters to an ODE, when using an adaptive step-size integrator.

OPTI includes a number of derivative free, global optimizers, including NOMAD (my favourite), PSwarm (parallelizable) and global solvers within NLopt.

## Example 1: Small GNLP
Lets revisit the non-convex problem from the NLP section:

![plot ex2nlpb](/img/opti/plot_ex2nlpb.png)

Previously we attempted to solve this problem using IPOPT, a local optimizer. Let's compare IPOPT's solution with those obtained from a selection of global optimization routines:

```matlab
% Objective
fun = @(x) norm([x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
                 x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))]);

% Bounds
lb = [-4;-4]; 
ub = [4;4];

% Initial Solution
x0 = [-4;-4];

% Build optiprob structure (intermediate structure)
prob = optiprob('fun',fun,'bounds',lb,ub,'x0',x0);

% Setup solver options
opts1 = optiset('solver','nomad');
opts2 = optiset('solver','pswarm');
opts3 = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
opts4 = optiset('solver','ipopt','warnings','off');

% Build OPTI objects
Opt1 = opti(prob,opts1); 
Opt2 = opti(prob,opts2); 
Opt3 = opti(prob,opts3); 
Opt4 = opti(prob,opts4); 
```

Note in this example we have used the intermediate function `optiprob` to create a common problem structure. The syntax for this function is the same as for `opti`, and can be more useful for larger problem descriptions.

Now let's solve each problem, and see what results we got:

```matlab
% Solve the GNLP problems
[x1,fval1] = solve(Opt1,x0);
[x2,fval2] = solve(Opt2,x0);
[x3,fval3] = solve(Opt3,x0);
[x4,fval4] = solve(Opt4,x0);
```

![plot ex1gnlp](/img/opti/plot_ex1gnlp.png)

The above plot shows each optimizer has in fact returned a different solution. Remember this problem has multiple *global* minima at fval = 0, thus solutions obtained by NOMAD and PSwarm appear valid (within the default tolerances). NLopt has found a reasonable solution, while IPOPT has fallen into a local minima, as expected.

You may obtain different results from PSwarm as it uses a randomizing pattern as part of it's algorithm. The other solvers are all deterministic.

## Example 2: Non-convex Polynomial {#scip-global}
The following example is a 6th order polynomial bounded between -1.5 and 1.5:

![ex gnlp poly](/img/opti/ex_gnlp_poly.png)

For this problem we will use a new feature of OPTI (starting from v1.73) which allows general MATLAB functions to be parsed into an algebraic description and then solved using SCIP. This allows SCIP to solve general NLPs and MINLPs to global optimality, provided an algebraic description of the system can be generated. The following code shows how SCIP can be substituted as the NLP solver and called just like any other NLP or GNLP solver:

```matlab
% Objective
fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;

% Bounds
lb = -1.5; 
ub = 1.5;

% Initial Solution
x0 = 0;

% Options
opts = optiset('solver','scip');
 
% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)
 
% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)

% Plot Solution within Problem Bounds
plot(Opt,[lb ub])
```

Using SCIP the above code will find the global minimum at `x` = -1.19, while local solvers such as IPOPT will typically get stuck in the local minimum at `x` = 0.486. Obviously moving `x0` will closer to the optimum will help local solvers, but this is not so easy to see in higher dimensions!

The new SCIP interface will parse any normal MATLAB function (or anonymous function) as long as the functions are limited `exp, log, log10, abs,` and `sqrt`. In addition, your function must be deterministic (no random or stochastic elements), and must not use any relational operators or trigonometric functions. Indexing and vectorizing up to 2D is supported. See the [White Box Optimization](../../guides/advanced/white-box.md) page for more examples.
