---
title: "Dynamic Optimization Settings"
slug: "/examples/dynamic-optimization/opts/"
---

This section will detail each of the options available via `optidynset`. These are options applicable only to solving dynamic optimization problems. To view all options, together with a summary, use:

```matlab
>> optidynset
```

Note all options are covered by relevant examples on the [DNLS Examples page](./dnls.md).

## integrator
The ODE integrator choice can make a HUGE difference when solving these problems, so ensure you know which integrator works best before trying to optimize the problem! `optidynset` allows the use of any of the MATLAB integrators, simply supply the name as a string (e.g. `'ode45'`, `'ode15s'`, etc). 

The default integrator is `ode45`.

## sensitivity
In order to obtain the gradient of the optimization problem with respect to the parameters, OPTI uses the sensitivity equation to generate the required derivatives. The sensitivity equation itself requires partial derivatives of the ODE function with respect to both states and parameters. 

The sensitivity option allows you to specify how you would like OPTI to approximate or use these derivatives. Options are as follows:

- `'ND'` - Numerical Differentiation, approximate using mklJac.
- `'AD'` - Automatic Differentiation, solve using autoJac.
- `'User'` - Use user supplied partial derivatives.
- `'None'` - Don't use sensitivity equation, optimization objective gradient is approximated by using mklJac across the integrator (typically very bad, but OK for gradient free solvers).
- Empty [] - Will use `'ND'` if a user supplied partial derivative function is not supplied, otherwise will use `'User'` if it is supplied.

The default is Empty.

## dfdz
Partial derivative function of the ODE with respect to the state variables (z). If not supplied it will be approximated using mklJac (ND).

## dfdp
Partial derivative function of the ODE with respect to the parameters (p). If not supplied it will be approximated using mklJac (ND).

## stateIndex
Indices (either numerical or logical) of the states to fit measured data to.

## initialT
Initial time to start the integrator. By default OPTI will use the first measurement time stamp, but this allows option allows you to solve from an earlier time.

## odeMaxTime
The integrator can get stuck on especially hard problems, so rather than wait, set a timeout on the integrator. The default is 30s, but you may like to considerably shorten this (unless your problem really does take 30s to solve!).

## odeOpts
Allows you to pass MATLAB integrator options (`odeset`) to your selected integrator. Do not pass a Jacobian or Mass Matrix, as these are not supported when solving with OPTI.
