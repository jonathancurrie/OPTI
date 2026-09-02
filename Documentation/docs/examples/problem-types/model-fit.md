---
title: "Model Fitting and Parameter Estimation"
slug: "/examples/problem-types/model-fit/"
---

This page introduces a new model fitting / parameter estimation object in OPTI v2.10, `optifit`. The purpose of this object is to create a common routine for solving model fitting problems, while also providing a library of common structures suitable for model identification.

The following examples introduce how to use `optifit` for both curve and surface fitting problems.

<small>**Note**</small> the `optifit` function does not currently solve dynamic system parameter regression problems. See the [DNLS](../dynamic-optimization/dnls.md) page for ideas on solving these problems with the `opti` object.

## The `optifit` Object
`optifit` is designed to be as simple as possible; just supply the predictor (or independent/experimental) variables as `xdata`, the outcome (or dependent) variables as `ydata` and the model as a string or function handle to fit, and `optifit` will attempt to solve the optimum parameters:

```matlab
% Curve fitting with optifit
>> ofit = optifit(xdata,ydata,model)
```

Returned from the above example is an `optifit` object, which contains the optimum parameters as the class property `theta`:

```matlab
% Access optimum model parameters
>> ofit.theta
```

If you are fitting a surface, simply concatenate `ydata` and `zdata` into a cell array as the second argument:

```matlab
% Surface fitting with optifit
>> ofit = optifit(xdata,{ydata,zdata},model)
```

The full calling form of `optifit` is:
```matlab
>> ofit = optifit(xdata,ydata,model,theta0,lb,ub,wts,opts)
```

where `theta0` is an initial parameter guess (required for fitting function handles), `lb` and `ub` are parameter bounds (optional), `wts` are data fitting weights (optional, the default is 1), and `opts` allows [`optiset`](../../guides/advanced/opts.md) options to be supplied to the object.

## Example 1: Curve Fitting A
To jump straight into it, here is an example solving a curve fitting problem:

```matlab
% Fitting Data
xdata = (1:40)';
ydata = [5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237,...
         2.7002,2.8781, 2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791,...
         1.6686, 1.6232, 1.571, 1.6057,1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129,...
         1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,0.96518, 1.2129, 1.2003, 1.0743]';

% Automatically identify model structure and solve optimal parameters
ofit = optifit(xdata,ydata,'auto')

% Plot the Result
plot(ofit)
```

When inspecting an `optifit` object, it will automatically print the solution statistics:
```c
---------------------------------------------------------------------------------
OPTI Fit Statistics
 R^2:             0.994748
 Adjusted R^2:    0.994464
 RMSE:            0.0983368
 SSE:             0.357795
 Model Structure: t1*exp(xd*t2) + t3 [exp2]
 Model Form:      Constant (With Intercept)
 Solver Status:   'Converged'

Nonlinear Least-Squares Analysis Of Variance:
 Source          DF     Sum of Squares      Mean Square       F Value       p Value
 Model            3            277.075          92.3582       3503.74   6.70698e-43
 Error           37           0.357795       0.00967013
 Uncorrected     40            277.432

Nonlinear Least-Squares Confidence Interval & Coefficient Statistics:
 Parameter      Estimate (   95% CI)     Std Error       t Value       p Value
  t1         5.37592 ±(    0.139714)     0.0689539        77.964   1.19963e-42
  t2      -0.0989081 ±(   0.0063353)     0.0031267      -31.6334   2.17271e-28
  t3         1.02149 ±(   0.0801745)      0.039569       25.8154   2.91147e-25
---------------------------------------------------------------------------------
```

Noting one of the most important details printed above is the model structure identified, which in this case is a three parameter exponential function: `t1*exp(xd*t2) + t3`, where `t` standards for `theta` and `xd` stands for `xdata`.

The overloaded `optifit` `plot` command will plot the data, the fitted solution and the confidence bounds (functional, simultaneous):

![ofit ex1](/img/opti/ofit_ex1.png)

## Example 2: Curve Fitting B
Here is another example where `optifit's` automatic structure identification is used:

```matlab
% Amount of Substrate
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; 
% Reaction Rate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; 

% Automatically identify model structure and solve optimal parameters
ofit = optifit(n,r,'auto')

% Plot the Result
plot(ofit)
```

This time a three parameter rational is chosen as the best structure. The resulting fit is shown below.

![ofit ex2](/img/opti/ofit_ex2.png)

## `optifit` Model Library
So far we have seen that `optifit` appears to work, both at solving reasonable parameters, but also in identifying reasonable model structures. The algorithm used in identifying the structure is a collection of heuristics, thus is certainly not robust. Therefore it may be that specifying a particular structure either provides better results, or based on the underlying physics, you know what the model structure should look like. The following table lists the current model structures built into `optifit`:

### Curve Fitting
<small></small>

|  |  |
| --- | --- |
| Name | Structure |
| `'poly1'` | ` t1*xd + t2` |
| `'poly2'` | ` t1*xd^2 + t2*xd + t3` |
| `'poly3'` | ` t1*xd^3 + t2*xd^2 + t3*xd + t4` |
| `'poly4'` | ` t1*xd^4 + t2*xd^3 + t3*xd^2 + t4*xd + t5` |
| `'poly5'` | ` t1*xd^5 + t2*xd^4 + t3*xd^3 + t4*xd^2 + t5*xd + t6` |
| `'power1'` | ` t1*xd^t2` |
| `'power2'` | ` t1*xd^t2 + t3` |
| `'exp1'` | ` t1*exp(t2*xd)` |
| `'exp2'` | ` t1*exp(t2*xd) + t3` |
| `'exp3'` | ` t1*exp(t2*xd) + t3*exp(t4*xd)` |
| `'rat01'` | ` t1 / (xd + t2)` |
| `'rat02'` | ` t1 / (xd^2 + t2*xd + t3)` |
| `'rat03'` | ` t1 / (xd^3 + t2*xd^2 + t3*xd + t4)` |
| `'rat11'` | ` (t1*xd + t2) / (xd + t3)` |
| `'rat12'` | ` (t1*xd + t2) / (xd^2 + t3*xd + t4)` |
| `'rat13'` | ` (t1*xd + t2) / (xd^3 + t3*xd^2 + t4*xd + t5)` |
| `'rat21'` | ` (t1*xd^2 + t2*xd + t3) / (xd + t4)` |
| `'rat22'` | ` (t1*xd^2 + t2*xd + t3) / (xd^2 + t4*xd + t5)` |
| `'rat23'` | ` (t1*xd^2 + t2*xd + t3) / (xd^3 + t4*xd^2 + t5*xd + t6)` |
| `'sin1'` | ` t1*sin(t2*x + t3) + t4` |
| `'sin2'` | ` t1*sin(t2*x + t3) + t4*sin(t5*xd + t6) + t7` |
| `'sin3'` | ` t1*sin(t2*x + t3) + t4*sin(t5*xd + t6) + t7*sin(t8*xd + t9) + t10` |
| `'auto'` | Will attempt to find the best structure from the above (up to 2nd order). |

### Surface Fitting
<small></small>

|  |  |
| --- | --- |
| Name | Structure |
| `'poly11'` | ` t1*xd + t2*yd + t3` |
| `'poly12'` | ` t1*yd^2 + t2*xd*yd + t3*xd + t4*yd + t5` |
| `'poly13'` | ` t1*yd^3 + t2*yd^2 + t3*xd*yd^2 + t4*xd*yd + t5*xd + t6*yd + t7` |
| `'poly21'` | ` t1*xd^2 + t2*xd*yd + t3*xd + t4*yd + t5` |
| `'poly22'` | ` t1*xd^2 + t2*yd^2 + t3*xd*yd + t4*xd + t5*yd + t6` |
| `'poly23'` | ` t1*yd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9` |
| `'poly31'` | ` t1*xd^3 + t2*xd^2 + t3*xd^2*yd + t4*xd*yd + t5*xd + t6*yd + t7` |
| `'poly32'` | ` t1*xd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9` |
| `'poly33'` | ` t1*xd^3 + t2*yd^3 + t3*xd^2 + t4*yd^2 + t5*xd^2*yd + t6*xd*yd^2 + t7*xd*yd + t8*xd + t9*yd + t10` |

In all models presented, `optifit` will automatically select a good initial parameter guess (using internal heuristics), as well as having the gradient of the model in its database. This enables it to have a fair chance of finding good parameters, if of course the model is a sensible choice for the data.

## Example 3: Fitting A Custom Nonlinear Function
If you do not want to use one of the built in models above, then you are free to pass a function handle of your own model. The function handle must be of the form `model(theta,xdata)`, or for surface fitting `model(theta,xdata,ydata)`. 

The example below is one of D. Himmelblau's data fitting examples:

```matlab
% Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure 
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate

% Model
Rxn = @(theta,xdata) theta(1)*xdata./(1+theta(2)*xdata);

% Fit Model Parameters
ofit = optifit(p,r,Rxn,[5e-3 2e-2])

% Plot the Result
plot(ofit)
```

The above code will automatically attempt to derive an analytical expression for the model gradient, build an OPTI object, select the best solver, and then solve the problem and return the solution as an `optifit` object. The plot below shows the model and 95% confidence bounds.

![ofit ex3](/img/opti/ofit_ex3.png)
