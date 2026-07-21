---
title: "Dynamic System Parameter Estimation"
slug: "/examples/dynamic-optimization/dynamic-system-parameter-estimation/"
---

Available from OPTI v2.00.

## Problem Overview
A dynamic parameter estimation problem aims to solve for the unknown parameters of a dynamic model, supplied as a system of Ordinary Differential Equations (ODEs). This problem can be setup as a standard Nonlinear Least Squares problem, however the nonlinear function to be fitted involves solving an ODE using a numerical integration scheme. OPTI calls this a Dynamic Nonlinear Least Squares (DNLS) problem.

The standard dynamic model form is shown below:

![def dnls](/img/opti/def_dnls.png)

`t`     - time  
`t0`    - initial time  
`z`     - state variable (noting `x` is reserved for optimization variables)  
`p`     - dynamic model parameter  
`z0`    - state initial condition  

Optionally, the problem can be posed to also (or alternatively) solve for unknown initial conditions. Therefore the solution vector (theta) is defined as follows:

![def dnls theta](/img/opti/def_dnls_theta.png)

The complication with DNLS problems is obtaining accurate gradients of the model with respect to the parameters, as the objective function also includes the ODE integrator. OPTI solves this problem by including the sensitivity differential equations as part of the solution process, which provide the gradient together with the objective.

### Importance of a good Initial Guess
These problems are incredibly sensitive to the initial solution guess (`theta0`) thus taking time (even manually if required) to 'tune' your initial guess is normally required. Even 1% away from the exact solution can result in substantial errors on tricky problems!

### Save Your Work!
Currently it is possible for the ODE integrator to get 'stuck' when solving these problems. If this happens, you may not be able to Ctrl-C to cancel execution. Therefore I suggest you save your work before solving these problems, as you may have to kill the MATLAB process in order to get 'unstuck'.

### Note on the Examples
The following examples demonstrate this beta functionality of OPTI, and are required to be completed in order. Each example will typically build on the functionality entered in the previous example.

## Example 1: Solving for Dynamic Model Parameters
Consider the following three state model with two unknown parameters:

![ex1 dnls](/img/opti/ex1_dnls.png)

To begin with, we will solve this model using 'known' parameters, to generate some fitting data:

```matlab
% ODE System
ode = @(t,z,p) [-p(1)*z(1) + 4; 
                2*z(1) - p(1)*z(2) + 5; 
                -4*z(1) - 2*z(2)*z(3) - p(2)];
 
% Initial Conditions
z0 = [-1.5;1.25;1];

% True Parameter Values
p = [2.345;1.1];
 
% Generate Fitting Data
tm  = 0:0.1:2;                 %measurement times
odeInt = @(t,z) ode(t,z,p);    %ODE function for ODE45
[~,zm] = ode45(odeInt,tm,z0);  %Solve ODEs
```

For this example we will have perfect data, but remember this is not often the case! To create this problem OPTI introduces two new arguments to the OPTI constructor, `ode` for supplying the ODE function, and `z0` for supplying the initial state vector. `data` functions the same as in NLS problems, while `theta0` is an alias for `x0`.

```matlab
% Starting Guess
theta0 = [1;0.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0)
```

Once the OPTI object has been created, the following information will be shown within the MATLAB command window:

```matlab
------------------------------------------------------
Dynamic Parameter Estimation Problem (DNLS)
 min sum[( int[F(t,z,p)] - ydata ).^2]
------------------------------------------------------
   Problem Properties: 
# Decision Variables:        2
  # Parameters:              2
  # Initial Conditions:      0
# Data Points:              63
------------------------------------------------------
  Solver Parameters:
Solver:                    MKLTRNLS
ODE Integrator:            ODE45
Sensitivity DFDZ:          Numerical Differentiation
Sensitivity DFDP:          Numerical Differentiation
------------------------------------------------------
```

OPTI has identified that our problem has two decision variables (two parameters to estimate), as well as 63 data points to fit the system parameters to (21 in each measurement vector x 3 states). OPTI has also defaulted to using the Intel Trust Region Nonlinear Least Squares solver, as this has been shown to perform well on these problems. In addition, MATLAB's `ode45` has been chosen as the default integrator.

The final two lines show the sensitivity equation partial derivatives are currently being calculated using numerical differentiation (via finite difference). We will see in a later example how to provide analytical derivatives for these terms, however finite difference typically works quite well.

To solve the problem it is as simple as calling `solve`:

```matlab
% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

Returned in the solution vector `theta` will be the same values as `p`, indicating the problem has been solved successfully. OPTI can also plot the solution, which is shown below. Smooth lines show the model response based on the solution `theta`, while grey circles show the measured data.

![plot ex1dnls](/img/opti/plot_ex1dnls.png)

## Example 2: Solving for Initial Conditions
For this example we will use the model and fitting data from Example 1, but also solve for initial conditions of states 2 and 3. Using OPTI, simply replace the 'unknown' initial conditions in `z0` with NaNs:

```matlab
% Indicate Initial Conditions to Solve for Using NaNs
z0 = [-1.5;NaN;NaN];
```

Then add the relevant initial condition guesses to `theta` and resolve the object:

```matlab
% Starting Guess [p;z0]
theta0 = [1;0.5;0.5;0.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

OPTI will plot grey squares to indicate the values of the solved initial conditions. As expected for this problem, they line up exactly over the measurement points.

![plot ex2dnls](/img/opti/plot_ex2dnls.png)

By estimating the initial conditions we are effectively allowing the ODE solution to *not* pass through the initial measurement point (assuming it is the same as the previously specified initial condition). This can result in a much better fit, at the sacrifice of not matching the initial point. In addition, there may be problems where you do not have data for all states at the initial point, thus estimating the initial condition(s) is required. Example 4 will demonstrate this feature.

## Example 3: Bounding Parameter Values
Quite often you will know the approximate range of the parameters and/or initial conditions you are fitting. OPTI can use this information as bounds in order to help the optimizer regress the parameters. For this example we will use the model and fitting data from Example 2, but also add bounds to the problem.

```matlab
% Add Bounds
lb = [1.5;0;0.1;0.1];
ub = [3.5;2;1.5;1.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'bounds',lb,ub,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)
```

Adding finite (and realistic) bounds on all variables will typically speed up the optimization process, and is generally recommended. In this example the solver uses one less function evaluation and takes one less iteration.

In addition, you can also specify linear inequality and equality constraints when solving a DNLS.

## Example 4: Fitting Selected States
Once again we will use the model and fitting data from Example 1, but this time we are only going to fit states 2 and 3. This is an example of when the system of differential equations contains one or more intermediate states, but we are actually only interested in fitting to output states where we have measured data. 

This examples introduces a new function, `optidynset`, which functions in a similar fashion to `optiset`. The difference is the `'dyn'`, indicating it is for setting dynamic optimization settings. `optidynset` allows the user to customize the integrator and derivatives among other settings.

```matlab
% Set Dynamic Options (specifying states of interest)
dopts = optidynset('stateIndex',[2 3]);

% Add to General OPTI Options
opts = optiset('display','iter','dynamicOpts',dopts);

% Index Measured Data
zm23 = zm(:,[2 3]);
```

Then as per standard OPTI problems, add the options structure to the OPTI argument list:

```matlab
% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm23,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

The plot below indicates OPTI has corrected fitted the parameters and initial conditions, even with only states 2 and 3. Remember this may not always be the case as you must always have enough information (i.e. non-zero terms in the derivatives of the ODE function) in order for OPTI to correctly identify theta.

![plot ex3dnls](/img/opti/plot_ex3dnls.png)

## Example 5: Multiple Sampling Rate Measurements
Sometimes you may have a set of fitting data (measurements) that have states (i.e. outputs) that have been sampled at different rates. This could be because temperature and flow rate is quick and easy to measure, while concentration or viscosity might take longer.

To illustrate, we will revisit the model from Example 1, but generate some more interesting fitting data.

```matlab
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm2  = 0:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
[~,zm2] = ode45(odeInt,tm2,z0); zm2 = zm2(:,2);
% Note 0 is required below just to match initial condition in this example
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3(2:end,3); %drop first point

% Group measurements and time stamps in cell arrays
tm_multi = {tm2;tm3};
zm_multi = {zm2;zm3};
```

When OPTI sees data specified as cell arrays (and is solving a DNLS), it will automatically setup the problem to correctly match the measured data against the ODE solution.

```matlab
% Return z0 to include NaNs of estimated initial conditions
z0 = [-1.5;NaN;NaN];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

The resulting plot shows that OPTI has once again correctly identified the parameters and initial conditions. State 3 (in red) passes through the less frequent measurement points, as well correctly solving for the initial condition (noting there was no measurement point at t = 0).

![plot ex4dnls](/img/opti/plot_ex4dnls.png)

## Example 6: 'Winding Back' the Integrator
Following on from the previous example, there may be a situation where you do not know any measurements at t = 0, but still wish to solve for initial conditions of the states at t=0. 

We will once again create some more interesting measurement data, but this time all our measurements will be at t > 0:

```matlab
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm1  = 0.5:0.1:2;              %state 1
tm2  = 1:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
tm_multi = {tm1;tm2;tm3};
zm_multi = {zm1;zm2;zm3};
```

As before, we have started each simulation above from t = 0 to recreate the 'real' data, however only included data at the designated measurement times. By default OPTI will start the ODE integrator at the first measurement time stamp, however we override it using `optidynset` and setting the 'initialT' setting to 0:

```matlab
% We will need to estimate all initial conditions in this problem
z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
theta0 = [1;0.5;0.5;0.5;0.5];

% Tell OPTI we want to solve from t = 0
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

The plot below shows the same solution as in the other examples (expected), but with the new 'interesting' data points. We also see OPTI has correctly identified all three initial conditions. 

![plot ex6tdnls](/img/opti/plot_ex6tdnls.png)

If you know one or more initial conditions at t = 0 you can enter them as per normal, but OPTI will provide a warning indicating that it is expecting these initial conditions are at t = 0. This is to alert users who supplying initial conditions of the first measurement point (rather than the start time of the ODE integrator).

## Example 7: Repeated Measurements
Good experimental design dictates that you should always have at least a few repeated measurement points. This ensures your measurement data is repeatable and realistic. 

We will once again create some more interesting measurement data, but this time have a few deliberate repeated data points:

```matlab
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% New ODE with modified parameters (representing a different run)
odeIntB = @(t,z) ode(t,z,p*0.85); 

% Measurement Times for each State
tm1  = 0.5:0.1:2;              %state 1
tm1b = 1:0.25:2;               %state 1 2nd run
tm2  = 1:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1]  = ode45(odeInt, [0 tm1], z0);  zm1  = zm1((2:end),1);
[~,zm1b] = ode45(odeIntB,[0 tm1b],z0);  zm1b = zm1b((2:end),1); %2nd run
[~,zm2]  = ode45(odeInt, [0 tm2], z0);  zm2  = zm2((2:end),2);
[~,zm3]  = ode45(odeInt, [0 tm3], z0);  zm3  = zm3((2:end),3);
[~,zm3b] = ode45(odeInt, [0 tm3], z0*0.75); zm3b = zm3b((2:end),3); %2nd run

% Group measurements and time stamps in cell arrays
tm_multi = {[tm1 tm1b];tm2;[tm3 tm3]}; %concatenate repeated points
zm_multi = {[zm1;zm1b];zm2;[zm3;zm3b]};
```

The code above simulates multiple runs of an experiment, collecting data at different rates for each state, as well as multiple runs for states 1 and 3. In the second run of state 1 the system parameters are slightly different (perhaps different ambient conditions), while in the second run of state 3 the initial conditions are slightly different (e.g. not started from the same point). Alternatively, the first run might be the 'incorrect' run (hard to ever know!). Each additional measurement vector is simply concatenated to the existing state measurements.

To solve this problem, the same code as the above example is used:

```matlab
% We will need to estimate all initial conditions in this problem
z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
theta0 = [1;0.5;0.5;0.5;0.5];

% Tell OPTI we want to solve from t = 0
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

As we now have extra data points not on the 'perfect' original trajectories, the solution is now quite different from previous examples. The plot below shows the optimizer has found a solution that appears to satisfy the least-squares criterion.

![plot ex7dnls](/img/opti/plot_ex7dnls.png)

## Example 8: Supplying Analytical Derivatives
Using a numerical approximation of the ODE partial derivatives provides a good starting point for solving these problems, but it can be slow, and on particularly numerically sensitive problems, result in errors. OPTI does provide an automatic differentiation option, but this is even slower than mklJac. Therefore substantial speed-ups can be achieved by supplying analytical derivative expressions.

The sensitivity differential equation used by OPTI takes the following form:

![def dnls S](/img/opti/def_dnls_S.png)

where

![def dnls J](/img/opti/def_dnls_J.png)

As you can see we require derivatives of the ODE with respect to both state variables and parameters. For small problems (such as Example 1), these can be solved easily using the MATLAB Symbolic Toolbox:

```matlab
% Declare Symbolic Variables
syms p1 p2 z1 z2 z3

% Symbolic ODE RHS
ODE = [-p1*z1 + 4; 
       2*z1 - p1*z2 + 5; 
       -4*z1 - 2*z2*z3 - p2];

% Solve Jacobians
dfdz_sym = jacobian(ODE,[z1 z2 z3])
dfdp_sym = jacobian(ODE,[p1 p2])
```

Alternatively, if your problem is entered in standard OPTI syntax (i.e. using the variables `t`, `z` and `p`, is *not* vectorized and you have the Symbolic Toolbox installed, you can use the following routine to automatically generate the required derivatives:

```matlab
% ODE System
 ode = @(t,z,p) [-p(1)*z(1) + 4; 
                 2*z(1) - p(1)*z(2) + 5; 
                 -4*z(1) - 2*z(2)*z(3) - p(2)];

% Simple Symbolic + String Parsing Derivative Generation
[dfdz,dfdp] = symDynJac(ode)
```

Either way, the derivative expressions can be entered as MATLAB functions for OPTI to use when solving the problem:

```matlab
% Analytical Derivative Expressions
dfdz = @(t,z,p) [-p(1) 0,       0;
                 2,    -p(1),   0;
                 -4,   -2*z(3), -2*z(2)];
dfdp = @(t,z,p) [-z(1), 0;
                 -z(2), 0;
                 0,    -1];

% Set Dynamic Options
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0);

% General OPTI Options
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);

% Create OPTI Object
 Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
            'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
```

On my computer I obtain a 3x speedup using the analytical derivative functions. Note we have also opted to enable the new OPTI derivative checker (via optiset), which will automatically validate the supplied derivatives.

## Example 9: Stiff Problems
So far our example problem is quite well behaved, and we haven't needed to use a special purpose integrator. Let's now look at the ODE equation for a flame, with an added parameter:

![ex9 dnls](/img/opti/ex9_dnls.png)

As with the first example, we define our ODE and generate some measurement data:

```matlab
% Flame ODE System
ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;
 
% Initial Condition (Controls Stiffness)
z0 = 0.001;

% True Parameter
p = 1.5;
 
% Generate Fitting Data
tm  = [0 65:70 80:0.1:90 100]*10;  %measurement times
odeInt = @(t,z) ode(t,z,p);        %ODE function for ODE15S
[~,zm] = ode15s(odeInt,tm,z0);     %Solve ODE
```

To demonstrate why the sensitivity equations <small>*can*</small> be useful (sometimes it is faster and more robust without!), we will solve the problem using the default integrator (`ode45`), with and without the sensitivity equations.

```matlab
% Starting Guess
theta0 = 1;

% Set Options (use ode15s)
dopts = optidynset('sensitivity','none');
opts1 = optiset('display','iter','iterfun',@optiplotlogfval,'dynamicOpts',dopts);
opts2 = optiset('display','iter','iterfun',@optiplotlogfval);

% Create OPTI Objects
OptNoDer = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts1)
OptWDer = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts2)

% Solve
[theta,fval,exitflag,info] = solve(OptNoDer)
[theta,fval,exitflag,info] = solve(OptWDer)

% Plot the Solution
subplot(121); plot(OptNoDer); title('No Sensitivity');
subplot(122); plot(OptWDer); title('With Sensitivity');
```

It is quite easy to see in the below plot that including sensitivity equations, while significantly slower to solve, means we are much closer to the true parameter.

![plot ex9dnls](/img/opti/plot_ex9dnls.png)

The reason this problem is especially 'delicate' relates to the fact it is stiff, meaning (roughly) that compared to the length of the solution time, one or more events happen very quickly (i.e. the flame igniting happens near instantaneously, compared to simply burning at steady state). If you zoom in on the final steady state in either plot you will see that the solution from `ode45` is bouncing around quite a bit. This bouncing, combined with the adaptive step-size functionality of the integrator means our derivatives by pure finite-difference (i.e. no sensitivity) are virtually worthless, hence the poor solution.

The solution to the bouncing (in this problem) is to use a stiff integrator, and a good one to try is `ode15s`. There is however a catch, stiff integrators use the Jacobian of system with respect to the states in order to choose an appropriate step size. By default, `ode15s` (and others) use finite-difference to determine the Jacobian (`dfdz` in our terminology). The problem now becomes apparent, if we use finite-difference (e.g. `mklJac`) to determine `dfdz` and thus include the sensitivity equations as part of the ODE, then MATLAB's finite-difference algorithm is attempting to find derivatives of our finite-difference derivatives. This results in a horrible mess, and can lead to exceptionally poor performance.

Therefore the rule is simple: if you use a stiff integrator you MUST supply analytical partial derivatives! If you don't, then turn off sensitivity (set to 'none'), and perhaps try a derivative free solver such as NOMAD.

Stiff Integrator with Analytical Derivatives:
```matlab
% Partial Derivatives
dfdz = @(t,z,p) 2*p(1)*z(1) - 3*z(1)^2;
dfdp = @(t,z,p) z(1)^2;

% Set Dynamic Options
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'integrator','ode15s',...
                   'sensitivity','User');

% General OPTI Options
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)
```

Stiff Integrator, Derivative Free Optimizer (in this case faster, not always the case)
```matlab
% Set Dynamic Options
dopts = optidynset('integrator','ode15s','sensitivity','None');

% General OPTI Options
opts = optiset('solver','nomad','display','iter','dynamicOpts',dopts);

% Create OPTI Object (always set finite bounds for derivative free)
Opt = opti('ode',ode,'data',tm,zm,'bounds',0.5,2.5,'z0',z0,...
           'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)
```
