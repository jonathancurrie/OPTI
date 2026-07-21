---
title: "What Is OPTI?"
slug: "/overview/"
---

OPTimization Interface (OPTI) Toolbox is a **free** MATLAB toolbox for constructing and solving linear, nonlinear, continuous and discrete optimization problems. A range of open source and academic solvers are supplied for the *Windows* user - no compilation required!

The toolbox is supported with demos, detailed documentation and implementation experience. The toolbox is released under the [BSD 3-Clause License](./project/license.md).

## What can I do with OPTI?
A brief outline of OPTI functionality is listed below. For more details, see the [Examples](./examples/index.md) section.

- Construct optimization problems in MATLAB and solve them using a range of supplied solvers
- Automatically detects the problem type being solved and uses the best solver available for the problem
- Use powerful open source solvers such as IPOPT, SCIP, NOMAD and others
- Perform solver benchmarking, plot optimization contours and validate the solution
- Read and write standard optimization files such as MPS and LP
- Read and solve AMPL problems

## A Typical User
A typical OPTI user will be:

- Comfortable using MATLAB (student level is fine)
- Understand the basics of optimization (how to define an objective and apply constraints)
- Working in Windows (sorry, no Linux or Mac support)
- Looking for a free optimization package

You may also be familiar with some of the following statements / questions which OPTI aims to solve:

- I'd like to use (insert solver here), but I'm not a strong C/C++ programmer
- I'm an OK programmer, but I just can't work out how to compile (insert solver here)
- Will a different solver give a better solution to my problem? 
- Can I try (insert solver here) without changing my problem?
- I'd like to solve Mixed Integer problems in MATLAB
- My problem is taking a long time to solve. Is there a faster solver which is still robust?

## Solvers
The main contribution of OPTI is providing the solvers compiled, with MEX interfaces. This allows a MATLAB user to call them just like any other MATLAB function.

All supplied solvers have been compiled from open source code obtained via the links listed in the [Solvers](./solvers/index.md) page. Where possible, all authors are aware of the inclusion in OPTI, and have OK'd this work. However if you use a solver in your research, please make sure you cite the original solver reference (and OPTI if you are feeling generous).

If you are using a solver in any commercial work please check the individual license restrictions of the solver.

## Further Reading

My [thesis](https://openrepository.aut.ac.nz/items/bb124ef8-830c-4446-80ec-5e5fc51f0bb0) also contains a huge amount of information on the development of OPTI , its motivation and some interesting case studies.
