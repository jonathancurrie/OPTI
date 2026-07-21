---
title: "Low Level Solver Interfaces"
slug: "/guides/advanced/low-level/"
---

It is entirely optional to use the OPTI class to use any of the solvers with this toolbox. For this reason most solvers are written with an easy to use interface if you choose to call them directly.

## MEX Interfaces & Wrappers
There are infact two interface levels for each solver, the lowest being the MEX interface itself, and the preferred level being a small MATLAB wrapper function over the top of the MEX interface.

To illustrate, consider the LP solver CLP:

### clp.mexw32 / clp.mexw64
- The MEX interface is a compiled C/C++ DLL which takes MATLAB arguments, passes them to the solver, runs it, then returns the result. Most error checking is performed at this level.

### opti_clp.m
- The MATLAB wrapper assigns solver specific defaults to unused arguments, checks the options against OPTI defaults, performs critical error checks, calls clp.mexw64, then creates the exitflag and information structure once the solver returns.

In order to further my suggestion to use the wrapper over the MEX interface itself (other than it may crash MATLAB or your computer if you get something wrong), lets consider the overhead introduced by the wrapper itself:

```matlab
% Load Large LP Test Problem
prob = coinRead('maros-r7.mps');

% Solve using OPTI Class
tic
solve(opti(prob,optiset('solver','clp')));
toc

% Solve using MATLAB Wrapper
tic
opti_clp([],prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub);
toc

% Solve using MEX Interface only
tic
clp([],prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub);
toc
```

On my PC I get the following times (run a number of times and averaged):

OPTI Class: 0.543s

MATLAB Wrapper: 0.522s

MEX Interface: 0.52s

Therefore for the sake of 2ms (for this problem) I recommend you use the wrapper at a minimum! And OPTI only uses < 50ms of overhead for this problem. Nonlinear problems will use substantially more overhead as OPTI runs a series of test function evaluations.

## Determining Calling Arguments
All MEX interface and wrapper functions are documented. For example:

```matlab
% MEX Interface Help
>> help clp

%% Wrapper Help
>> help opti_clp
```

Simply enter the solver name (as returned by optiSolver) as the argument to help to determine the calling arguments.

## Example 1: NLS
This example will show how to call the Intel Math Kernel Library (MKL) trust region nonlinear least squares solver using the low level interface. We will use the Rosenbrock banana function again, as this was originally a NLS problem.

```matlab
% Fitting Function
fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];

% Fitting Data
ydata = [0;0];

% Starting Guess
x0 = [-1.2;1];

% Call Low Level Interface
[x,fval] = opti_mkltrnls(fun,[],x0,ydata)
```

Note for many NLS solvers they have a built in finite difference algorithm so a gradient is not required. By default however OPTI will use it's own finite difference algorithm when a NLS solver is called via the class.

## Example 2: NLP {#ipopt}
Unfortunately IPOPT, NLOPT and BONMIN all require too many arguments to make a simple interface like those described so far. Therefore these use the MEX interface as the preferred calling method for low level use.

Consider the following NLP to solve with IPOPT:

![ex2 ll nlp](/img/opti/ex2_ll_nlp.png)

```matlab
% Objective & Gradient
fun = @(x) log(1+x(1)^2) - x(2);
grad = @(x)[(2*x(1))/(x(1)^2+1), -1];

% Nonlinear Constraint, Jacobian & Structure
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2 - 4;
cl = 0; cu = 0;
nljac = @(x) sparse([4*x(1)*(x(1)^2+1),2*x(2)]);
jacstr = @() sparse([1 1]);        

% Starting Guess
x0 = [2;2];

% Build Function Structure
funcs.objective = fun;
funcs.gradient = grad;
funcs.constraints = nlcon;
funcs.jacobian = nljac;
funcs.jacobianstructure = jacstr;

% Build Options Structure
opts.lb = -1*Inf(2,1);
opts.ub = Inf(2,1);
opts.cl = cl;
opts.cu = cu;
opts.ipopt.hessian_approximation = 'limited-memory';

% Call IPOPT
[x,output] = ipopt(x0,funcs,opts)
```

As you can see a bit more work is required to generate the required structures suitable for calling IPOPT. Using the OPTI class this is taken care of automatically, including approximating derivatives if they are not supplied. For comparison, here is the same problem using the OPTI class and automatically approximated derivatives:

```matlab
% Build OPTI Object
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'x0',x0)

% Solve
[x,fval] = solve(Opt)
```

For toy problems approximating derivatives is fine. However when solving real optimization problems I highly recommend exact derivatives. See the [1st Derivatives](./deriv1.md) section for more details on how to do this with OPTI.

## Summary
While calling the solvers without OPTI is fine, using the class does have some advantages. These include more detailed error checking, automatic conversions between solvers, a common calling syntax and automatic option and derivative settings.
