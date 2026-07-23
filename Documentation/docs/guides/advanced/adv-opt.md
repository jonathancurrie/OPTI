---
title: "Advanced Options"
slug: "/guides/advanced/adv-opt/"
---

This section will detail each of the options available via `optiset`. These are options which are generally applicable to most solvers, so are grouped for easy use. To view all options, together with a summary, use:

```matlab
>> optiset
```

Remember not all solvers support all options. To view which options are applicable to your selected solver, use:

```matlab
>> optiSolver('config')
```

## solver
You can override OPTI's solver choice by selecting it here. Enter a string containing the name of the solver to use. In some instances OPTI may still override your choice if your solver is not available or not compatible with the problem you have entered.

The default is `auto`, where OPTI will choose a solver based on your problem.

## display
When you start solving real problems you will want to monitor the solver's progress, and this option is the best way to do this! Using `display` you can enable the solver to print it's progress to the MATLAB command window. Select from `'off'` (no display), `'final'` (final summary of solver run) or `'iter'` (print iteration by iteration results).

The default is `off`.

## maxiter
Limit the number of iterations the solver can run. Note one iteration may contain multiple function, gradient and constraint calls. 

The default is 1500 iterations.

## maxfeval
Limit the number of objective function evaluations (nonlinear problems only) the solver can call. 

The default is 10000 evaluations.

## maxtime
Limit the maximum execution time of a solver. This is a very useful option as it allows you to return the best solution given the maximum time you are happy to wait. If the solver finishes before this time, then good!

The default is 1000 seconds (16.66 minutes).

## maxnodes
Limit the maximum nodes a mixed integer solver can explore. Each node typically represents at least one relaxed (continuous) problem that will be solved. 

The default is 10000 nodes.

## tolrfun
The desired relative convergence tolerance of the solver. This is defined as:

$$
\frac{\left| f_i - f_{i-1} \right|}{\left| f_i \right|} \leq \mathit{tol}_{\mathit{rel}}
$$

Where *f<sub>i</sub>* is the current objective value, and *f<sub>i-1</sub>* is the previous iteration's objective value. The default is 1e-7.

## tolafun
The desired absolute convergence tolerance of the solver. This is defined as:

$$
\left| f_i - f_{i-1} \right| \leq \mathit{tol}_{\mathit{abs}}
$$

Where *f<sub>i</sub>* is the current objective value, and *f<sub>i-1</sub>* is the previous iteration's objective value. The default is 1e-7.

## tolint
The absolute tolerance used to define whether a solution is an integer value. This is defined as:

$$
\left| x_r - x_z \right| \leq \mathit{tol}_{\mathit{int}}
$$

Where *x<sub>r</sub>* is the real valued solution, and *x<sub>z</sub>* is the integer value. The default is 1e-5.

## warnings
OPTI will always indicate via a warning whether it is making an assumption or converting the problem from your original form. From OPTI v1.80 there are now three levels of warnings, allowing you to specify the level detail you are interested in. Options are `'all'` (previously `'on'`), `'critical'` or `'none'` (previously `'off'`). Previous options are still compatible with newer versions of OPTI.

The default is `'critical'`. If you would like to try and speed up the optimization process, try changing the warning level to `'all'`, then address each of the warnings OPTI prints.

## iterfun
Most nonlinear solvers allow a callback method to be automatically called at every iteration (or function evaluation, depending on the solver). The callback function must have the following form:

```matlab
function stop = myCallback(iter,fval,x) 
```

Where the input arguments are: 

- `iter` - The current iteration or function evaluation
- `fval` - The current objective function value
- `x` - The current solution vector (not available from all solvers) 

And the output arguments are: 

- `stop` - Set this to true to stop the solver when the function returns. The default is false (continue). Do not attempt to set it to 0 or 1, only logical values are accepted (not available from all solvers). 

The callback function can be used to set customized stopping criteria, customized iteration display, or for plotting. OPTI provides a default plotting function for plotting the objective value vs iteration. This can be set as follows: 

```matlab
opts = optiset('iterfun',@optiplotfval)
```

## solverOpts
This field allows you to supply solver specific options via a structure. For example if you are using IPOPT and wish to enable the derivative test option, it can be passed via this field. The same concept applies for many other solvers such as MATLAB's `fmincon`, PSwarm, NOMAD, BONMIN and others.

To illustrate, consider we want to tell NOMAD to stop once a certain objective value is reached. This option is set via `nomadset` (<small>**note**</small> the naming strategy, it is always the solver name appended with 'set') via the option `f_target`:

```matlab
% Set NOMAD specific options
nopts = nomadset('f_target',1e-6);

% Group with all OPTI options
opts = optiset('display','iter','solverOpts',nopts)
```

Note there may be common options between the solver and optiset. In this case optiset options always take priority.

To see which solvers contain extra option methods, check the `sOpts` column when using `optiSolver('config')`.

## Overriding OPTI's Problem Identification {#probtype}
For 1D problems OPTI may not recognise the type of problem you have posed correctly. If this is the case, you can override OPTI's problem identification by specifying an extra parameter to `opti` or `optiprob`, `'probtype'`.

An example of when OPTI will get it wrong is when finding the root of a single variable function (SNLE). OPTI will incorrectly identify this as a Unconstrained Nonlinear Optimization (UNO) problem:

```matlab
% Equation to solve for zero
fun = @(x) x.^3-2*x-5

% Starting Guess
x0 = 2;

% Build OPTI Object, specifying problem type
Opt = opti('fun',fun,'x0',x0,'probtype','SNLE')

% Solve SNLE problem
[x,fval,exitflag,info] = solve(Opt)
```

This problem occurs because OPTI expects a SNLE to contain more than one equation. While the above problem is fine, it contains no distinguishing features from a UNO, thus you must manually specify the problem type in these instances.
