---
title: "Optimization Settings"
slug: "/guides/advanced/opts/"
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

## tolrfun {#tolrfun}
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

## iterfun {#iterfun}
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

## solverOpts {#s-opts}
This field allows you to supply solver specific options via a structure. For example if you are using IPOPT and wish to enable the derivative test option, it can be passed via this field. The same concept applies for many other solvers such as MATLAB's `fmincon`, PSwarm, NOMAD, BONMIN and others.

To illustrate, consider we want to tell NOMAD to stop once a certain objective value is reached. This option is set via `nomadset` (**note** the naming strategy, it is always the solver name appended with 'set') via the option `f_target`:

```matlab
% Set NOMAD specific options
nopts = nomadset('f_target',1e-6);

% Group with all OPTI options
opts = optiset('display','iter','solverOpts',nopts)
```

Note there may be common options between the solver and optiset. In this case optiset options always take priority.

To see which solvers contain extra option methods, check the `sOpts` column when using `optiSolver('config')`.

## dynamicOpts
This field allows you to supply settings specific to dynamic optimization via structure. This structure is generated using [`optidynset`](../../reference/options/optidynset.md) and is passed to `optiset` in the same way as solver options, as above.

## derivCheck {#deriv-check}
When supplying derivative functions to an optimizer it is essential they are correct. From v2.00 OPTI can now automatically verify your derivatives are correct by comparing the supplied function(s) (and sparsity patterns, if supplied) to an internal numerical approximation. OPTI will check both first and second derivatives, if supplied. If there are errors, OPTI will report which element it considers to have an error in the command window.

This option is enabled by setting it to `'on'`. By default it is `'off'`.
