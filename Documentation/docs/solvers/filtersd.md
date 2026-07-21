---
title: "FILTERSD"
slug: "/solvers/filtersd/"
---

## FilterSD [Supplied]
<small></small>
FilterSD solves smooth, twice differentiable, nonlinear programs. In addition, it contains both dense and sparse solvers, which allows for large-scale NLPs to be solved (provided you can supply a sparse Jacobian). FilterSD is a local solver only, but is pretty fast!

|  |  |
| --- | --- |
| Category | Open Source |
| Author | [Prof. Roger Fletcher](http://www.maths.dundee.ac.uk/fletcher/) |
| License | [Eclipse Public License](http://www.opensource.org/licenses/eclipse-1.0) |
| Home Page | [FilterSD Home Page](https://projects.coin-or.org/filterSD) |
| Download Page | [FilterSD Tarball Directory](http://www.coin-or.org/download/source/filterSD/) |
| MEX Interface | OPTI Version Supplied |
| Pre-Requisites | [MKL](http://software.intel.com/en-us/articles/intel-mkl/) |
| Version Supplied | 1.0 |

<small>**Note**</small> I have updated the interface to dynamically allocate memory based on a problem size estimate. However when using the sparse version for really big problems (10,000 variables +) I have seen some crashes. Please report this to me if you find the same issue!
