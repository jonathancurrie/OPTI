---
title: "DJACOBI"
slug: "/reference/utilities/djacobi/"
---

## Intel MKL DJACOBI [Supplied]
<small></small>
DJACOBI is the default numerical difference routine for all derivative based nonlinear optimizers in OPTI. It uses a centered-difference approach to ensure derivatives are approximated as accurately as possible.

|  |  |
| --- | --- |
| Category | Commercial |
| Manager | Intel |
| License |  |
| Home Page | [Intel MKL Home Page](http://software.intel.com/en-us/articles/intel-mkl/) |
| Download Page | [Intel MKL Download Page](http://software.intel.com/en-us/articles/intel-mkl/) |
| MEX Interface | OPTI Version Supplied |
| Pre-Requisites | [MKL](http://software.intel.com/en-us/articles/intel-mkl/) |
| Version Supplied | 11.3 Release 3 |

### Sequential Build
Note this function is built against the sequential MKL. This is due to the MEX interface callback to MATLAB not being thread safe. If anyone can work out how to make it thread safe, I can enable parallel finite-difference approximations!

### Intel MKL
The Intel Math Kernel Library (MKL) is a suite of powerful mathematical routines encompassing BLAS, LAPACK, Sparse BLAS, nonlinear solvers, and many other useful routines. It is one of the few packages I am really impressed with, and use it extensively in my own research. Academic pricing is really competitive, and I thoroughly recommend it!
