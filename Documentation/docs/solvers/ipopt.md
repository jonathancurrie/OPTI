---
title: "IPOPT"
slug: "/solvers/ipopt/"
---

## Interior Point Optimizer [Supplied]
<small></small>
IPOPT solves smooth, twice differentiable, nonlinear programs. While the objective need not be convex, IPOPT will only find local solutions. However for convex NLPs it is probably the best open source solver available!

|  |  |
| --- | --- |
| Category | Open Source |
| Manager | [Dr. Andreas Wächter](https://www.mccormick.northwestern.edu/research-faculty/directory/affiliated/waechter-andreas.html) |
| License | [Eclipse Public License](https://opensource.org/license/EPL-1.0) |
| Home Page | [IPOPT Home Page](https://projects.coin-or.org/Ipopt) |
| Download Page | [IPOPT Tarball Directory](http://www.coin-or.org/download/source/Ipopt/) |
| MEX Interface | Supplied with IPOPT (by  [Dr. Peter Carbonetto](https://profiles.uchicago.edu/profiles/display/37244)) |
| Pre-Requisites | [MA57'^1^'](https://www.hsl.rl.ac.uk/catalogue/ma57.html), [MUMPS'^1^'](./mumps.md),  [MKL'^2^'](http://software.intel.com/en-us/articles/intel-mkl/) |
| Version Supplied | 3.12.7 |

### Sparse Linear Solver
'^1^'As of OPTI Toolbox v1.71 the HSL solver MA57 has been dynamically linked to the IPOPT Interface. This has been achieved as MA57 is supplied with MATLAB, thus OPTI simply uses the version of MA57 already on your computer! MA57 appears to solve all problems via IPOPT faster than MUMPS, and should also be more robust.

'^2^'In addition, from OPTI v2.05 MKL PARDISO is available as a sparse linear solver. PARDISO is parallelized and can substantially speedup solving large problems. It is however a little 'sensitive', so ensure you have specified derivatives and sparsity patterns correctly. As below, select PARDISO via the linear_solver option.

<small>**Note**</small> it is possible to still use MUMPS. This is achieved by changing the 'linear_solver' option in ipoptset().

### Citation
A. Wächter and L. T. Biegler, "On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming," *Mathematical Programming* 106(1), pp. 25-57, 2006
