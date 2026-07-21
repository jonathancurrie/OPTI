---
title: "Dynamic System Parameter Estimation"
slug: "/examples/dynamic-optimization/dynamic-parameter-estimation/"
---

## Problem Overview
A dynamic parameter estimation problem aims to solve for one or more parameters of a dynamic model, supplied as a system of Ordinary Differential Equations (ODEs). This problem can be setup as a normal Nonlinear Least Squares problem, however the nonlinear function to be fitted involves solving the ODE, using a numerical integration scheme. OPTI calls this a Dynamic Nonlinear Least Squares (DNLS) problem.

The complication with DNLS problems is obtaining accurate gradients of the model with respect to the parameters, as the objective function also includes the ODE integrator. OPTI solves this problem using the established method of including the sensitivity differential equations as part of the solution process, which provide the gradient together with the objective.
