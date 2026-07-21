---
title: "White Box Optimization"
slug: "/guides/advanced/white-box/"
---

Nonlinear optimization using the solvers described in the [NLP](../../examples/problem-types/nlp.md) and [GNLP](../../examples/problem-types/gnlp.md) sections, is deemed *black-box* because as far as the optimizer is concerned, the internals of the nonlinear objective and constraint functions are unable to be exploited. This blackbox simply accepts a decision variable vector and returns an objective value or vector of constraint evaluations. The underlying equations, structure, and relationships between the decision variables (inputs) and objective or constraints (outputs) are effectively unknown, with the exception of the derivatives, if available, and even derivative information is limited.

In contrast, linear and quadratic optimizers have a rigid mathematical format which the solver is tailored to exploit and solve. Furthermore, because a linear or quadratic program is described natively using numerical matrices and vectors, the solver has a full problem description which includes all problem data. This means that higher levels of pre-processing, such as scaling, redundant constraint removal and bound tightening can be applied. In addition, the standard problem definition of a linear or quadratic program requires the problem to be convex, further assisting in tailoring the internals of the solver to exploit matrix properties (such as positive-definite systems of linear equations). This all combines to make linear and quadratic solvers typically much faster and more robust than their general nonlinear counterparts, and they are therefore capable of solving much larger problems within a reasonable time frame.

To be able to leverage the same level of problem information for a nonlinear problem requires an algebraic description of the problem supplied to the solver. This effectively *opens-up* the back-box so that the optimizer can see inside and thus results in the term *white-box*. By exploiting the algebraic description, a white-box solver is able to much more effectively pre-process the problem, which can vastly reduce the search space required. In addition, a white-box solver may recognise mathematical features of the supplied functions, such as identifying monomials or polynomials, linear, bilinear and multilinear relationships or other common nonlinear expressions such as log and exp. This enables the white-box solver to be able to exploit problem relaxations and solution heuristics, such as outer approximations, convex/concave envelopes, a Generalized Benders Decomposition (GBD) and/or search space cutting planes. In addition, due to the solvers access to the algebraic description of the problem, it can generate derivatives internally, thus there is no need to supply 1st or 2nd derivatives to a white-box solver.

The most powerful feature of a white-box solver though is by using a suite of mathematical rules based on the structure of the problem, it can **prove deterministically a global optimum to general nonlinear problems, including non-convex and integer problems**, provided that the problem meets certain formulation requirements.

OPTI provides interfaces to two white-box optimization solvers via MATLAB: the open-source (but only available free for academics) [SCIP](../../solvers/scip.md), and the commercial [BARON](../../solvers/baron.md).

## Reproducibility using White Box Solvers
Be aware that different versions of MATLAB, including different architectures (x86 vs amd64) can and will cause minor differences in the solutions found by these solvers. In very rare cases, they may even report a local solution as global, but I have only seen this once. You are more likely to simply see different numbers of iterations between different PCs and MATLAB versions, but the same solution at the end.

Also if your problem is particularly poorly scaled, these solvers will have a hard time. They are sensitive to numerical round-off, and double precision just doesn't always cut it.

## Problem Restrictions
The biggest hurdle to using a white-box solver is supplying an algebraic description of the optimization problem, and the restrictions this places on how your model can be formulated. The list below summarizes the main restrictions when using an OPTI interface to a white-box solver.

- *Code* - Both OPTI interfaces require that your *entire* optimization problem is written in MATLAB code. This means you cannot call MEX functions, Simulink models or other external routines.
- *Functions* - Only a small subset of nonlinear functions are permitted by the solvers. Generally these are `power`, `log`, `log10`, `exp`, and in some cases, `abs` and `sign`. Note all trigonometric functions, statistical functions, interpolation, `max`, `min` and other common functions are not supported.
- *Determinism* - Your optimization problem must be deterministic in order to be processed by a white-box solver, as well as to be solved to global optimality. This means you cannot use any conditional statements (`<`,`>`,`<=`,`>=`,`!=`,`==`) with decision variables, as these provide an alternate code path for the optimization problem that is not supported. In addition, you must not add any random numbers to your model, as during optimization they will be treated as a constant (the call to your objective/constraints is made once during parsing, hence only one set of random numbers will be used).

While the above list may seem overly restrictive, many of the examples in this documentation can still be solved, as will many industrial problems. In addition, you are still free to use `for` loops (provided they are deterministic), multiple `.m` files, sparsity and matrix operations (although best handled using BARON currently).

Furthermore, with a well posed (bounded) optimization problem, you are free to use polynomials and other supported functions to approximate nonlinear functions over the region of interest, allowing you to still solve your problems of interest.

## SCIP
[SCIP](../../solvers/scip.md) is an open-source optimization solver designed for solving continuous and integer problems with linear and quadratic constraints, thus problems up to [MIQCQP](../../examples/problem-types/miqcqp.md) can be reliably solved using SCIP. It also includes beta functionality for solving general nonlinear problems, and in my experience, for small-medium nonlinear problems it is quite fast and robust. For academics using the Academic version of OPTI, SCIP is included for free. Otherwise, you will need to take out a license. See the [SCIP](../../solvers/scip.md) solver page for information on licensing.

The OPTI interface to the nonlinear functionality in SCIP is quite rudimentary, with only a few MATLAB methods included, as well as being reasonably inefficient. However the interface is robust and self-validating and I have not had any complaints (so far). Use the following command to see what functions are available with SCIP:

```matlab
>> methods(scipvar)
```

Note that only matrix multiplication is supported, other matrix routines including matrix division are not supported. Sparsity is also ignored, so be careful with large problems. Finally the SCIP interface is not designed to allow you to concatenate numerical values with your decision variable vector. For fixes to most of the above, use BARON (discussed below).

### SCIP Algebraic Description
I designed the SCIP algebraic description using what I thought was a novel scheme I termed an 'instruction list'. It wasn't until my Ph.D. supervisor saw what I was doing that he said I had reinvented  [Reverse Polish Notation](http://en.wikipedia.org/wiki/Reverse_Polish_notation). Regardless, the code was written, and the interface worked satisfactorily. 

Two components work together to allow SCIP to understand MATLAB code; a MATLAB class, `scipvar`, and MEX function, `scip`. The MATLAB class takes your MATLAB function and converts it to (effectively) RPN, or as I called it, an instruction list, as shown in the below example:

```matlab
% SCIP variable vector
x = scipvar(2,1);
% Test Function
fun = @(x) log(1 + x(1)^2) - x(2);

%Evaluate Function using SCIP variables
fun(x)
```

The result is an instruction list:
```matlab
VAR   : 0
NUM   : 2
POW   : NaN
NUM   : 1
ADD   : 1
LOG   : NaN
VAR   : 1
SUB   : NaN
```

Where the left hand column contains the instructions, and the right hand column the data. Variables are indexed from 0, and NaNs are ignored as placeholders for functions. Reading from top to bottom, you should see the same operation as the test function.

This instruction list is then passed as an array to the SCIP MEX function, which then builds an "Expression Tree" within C++, which is used by SCIP and CppAD for solving the problem and generating derivatives.

### SCIP Example 1
The following problem is number 20 in the Hock-Schittkowski collection, and is the only problem IPOPT currently gets wrong in the OPTI NLP benchmark:

![ex1 scip](/img/opti/ex1_scip.png)

```matlab
% Objective
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) [x(1) + x(2)^2;                    %cl <= nlcon(x) <= cu
              x(1)^2 + x(2);
              x(1)^2 + x(2)^2 - 1];
cl = [0;0;0];
cu = [Inf;Inf;Inf];

% Bounds
lb = [-0.5;-inf];                               %lb <= x <= ub
ub = [0.5;inf];

% Create OPTI Object
opts = optiset('solver','scip','display','iter');
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'opts',opts)

% Solve the NLP problem
x0 = [-2;1] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)
```

Examining the pre-solving output of SCIP, we see:
```matlab
presolved problem has 4 variables (0 bin, 0 int, 0 impl, 4 cont) and 5 constraints
      5 constraints of type <quadratic>
```

Which shows, even though the objective contained a quartic, SCIP has managed to convert the problem to a problem with only quadratic constraints, thus it is well set up to solve the optimization problem. This is typical of SCIP, as the pre-solving step is particularly efficient.

Examining the solution output:
```matlab
SCIP Status        : problem is solved [optimal solution found]
Solving Time (sec) : 0.26
Solving Nodes      : 1
Primal Bound       : +3.81987165140230e+001 (1 solutions)
Dual Bound         : +3.81987165140230e+001
Gap                : 0.00 
  [nonlinear] <NonlinearObj0>: (( +100 power((<xvar1> - sqr(<xvar0>)), 2)) + power((1 -1 
                                  <xvar0>), 2))-1<nlobj>[C]  == 0;
violation: right hand side is violated by 6.20194767719795e-006 (scaled: 6.20194767719795e-006)
best solution is not feasible in original problem
```

We see the problem is solved (and quite quickly), with the Gap between the primal and dual problems 0.0, indicating a global solution has been found. However we also see a warning indicating the 'solution is not feasible in the original problem'. This can be due to a couple of reasons:

- Your objective function and/or constraints are not scaled appropriately, and numerical errors have caused the transformed, SCIP problem, to differ from the original user problem. In this case, the objective could do with a little rescaling to remove the violation (although it is tiny):
```matlab
fun = @(x) 10*(x(2)-x(1)^2)^2 + 0.1*(1-x(1))^2;
```
- We have both unbounded decision variables and constraints within our optimization problem. These can cause both problems converging (much slower and perhaps not at all), and numerical problems as automatic scaling may not be able to be applied. *Always* apply finite bounds to all variables and constraints, if possible.

### SCIP Example 2 {#scipopts}
Rather than use SCIP to find the global solution, it may be advantageous to leverage the automatic differentiation within SCIP to simply find a valid local solution, and therefore avoid the typically long solve times of the full global solution, and the complexity of entering derivatives for other solvers.

![ex1 nlp](/img/opti/ex1_nlp.png)

```matlab
% Objective
fun = @(x) log(1 + x(1)^2) - x(2);    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2; %Nonlinear Equality Constraint
cl = 4;
cu = 4;

% Create SCIP Options [limit to 1 solution]
sopts = scipset('scipopts',{'limits/solutions',1});

% Create OPTI Object
opts = optiset('solver','scip','display','iter','solverOpts',sopts);
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'ndec',2,'opts',opts)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)
```

In the above example `scipset` has been used to customize the options for SCIP, which are then passed to `optiset` and then to `opti`, as per normal OPTI usage. SCIP stops once the first feasible local solution is found, and returns the associated decision variable vector.

For a full list of SCIP options, consult the  [SCIP options list](http://scip.zib.de/doc/html/PARAMETERS.php), noting they are entered a little differently from most other solvers.

### SCIP Example 3
Note SCIP has no problem if you want to add integer constraints as well:

![ex3 scip](/img/opti/ex3_scip.png)

```matlab
% Objective
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

% Linear Constraints
A = [-1 1; 1 1];                                %rl <= A*x <= ru
rl = [-inf;5];
ru = [-1;5];

% Bounds
lb = [0;0];                                     %lb <= x <= ub
ub = [4;4];

% Integer Constraints
xtype = 'IC';

% Create OPTI Object
opts = optiset('solver','scip','display','iter');
Opt = opti('fun',fun,'lin',A,rl,ru,'bounds',lb,ub,'xtype',xtype,'opts',opts)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)
```

Running the above code SCIP solves just as quickly, and returns the globally optimal integer solution to the optimization problem.

## BARON
[BARON](../../solvers/baron.md) is a commercial optimization solver that is designed for solving large, non-convex, integer, nonlinear optimization problems, and is one the fastest (if not the fastest) solvers available for these problems. Together with the developer of BARON, I developed a MATLAB interface for BARON which could exploit the full functionality of BARON, as well as attempting to provide the most compatible interface for traditional MATLAB users. BARON is definitely worth trying if you are serious about large-scale nonlinear optimization.

BARON is available as a free demo for problems with up to 10 variables, and the interface is also freely available. Have a look on the [BARON](../../solvers/baron.md) page for more information on obtaining BARON. Note the BARON interface is directly compatible with OPTI, thus BARON can be easily selected as the solver.

To see what functions are available with BARON, use:
```matlab
>> methods(barvec)
```

Note while many more matrix operations are compatible with BARON, including sparsity and simple transformations, matrix division is not supported. BARON also does not support any extra nonlinear functions. The interface is still in active development, so if it doesn't do what you need, then send me an email and I will determine if it can be added.

Also be aware that BARON is not called as a MEX file, instead it is an executable which is passed the optimization problem as a text file, solves the problem, then produces a text solution file which is processed by the interface to return the results to MATLAB.

### BARON Algebraic Description
BARON does not actually convert your MATLAB problem to an algebraic description, rather it parses it to a text format that is compatible with the BARON parser. Following the same example as with SCIP: 

```matlab
% BARON variable vector
x = barvec(2,1);
% Test Function
fun = @(x) log(1 + x(1)^2) - x(2);

%Evaluate Function using BARON variables
fun(x)
```

The result is a string of the function:
```matlab
Eq : log(1 + x1^2) - x2
```

Noting the output is a string with the variables renamed, and appropriate transformations made if required. The string is then written to a BARON compatible problem file, the executable called from MATLAB, and the problem solved and solutions retrieved.

### BARON Example 1
Using the same example from SCIP Example 1 above:

```matlab
% Objective
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) [x(1) + x(2)^2;                    %cl <= nlcon(x) <= cu
              x(1)^2 + x(2);
              x(1)^2 + x(2)^2 - 1];
cl = [0;0;0];
cu = [Inf;Inf;Inf];

% Bounds
lb = [-0.5;-inf];                               %lb <= x <= ub
ub = [0.5;inf];

% Create OPTI Object
opts = optiset('solver','baron','display','iter');
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'opts',opts)

% Solve the NLP problem
x0 = [-2;1] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)
```

The problem is solved easily, and without any numerical warnings as reported by SCIP. Note the only difference required was the change of solver name, OPTI handled the conversion to the BARON interface succinctly. 

BARON can be used to solve all other problems presented on this page. These are left for you to try!

## Summary
As shown, SCIP and BARON are unique solvers which exploit the structure of an optimization problem, in a similar way LP and QP solvers do, but for nonlinear problems as well. This allows them to not only find better solutions faster and more robustly, but also prove a solution is a global solution, which is quite a feat!
