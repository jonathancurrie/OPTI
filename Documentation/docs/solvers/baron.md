---
title: "BARON"
slug: "/solvers/baron/"
---

## BARON [Interfaced]
<small></small>
BARON is a commercial global optimization solver for linear and quadratic and nonlinear programs, including those with mixed integer constraints. It deterministically proves a global optimum to all compatible problems, including non-convex and mixed integer varieties. BARON is the world-leader in global optimization technology and its performance will not disappoint! 

|  |  |
| --- | --- |
| Category | Commercial |
| Manager | [Prof. Nick Sahinidis](http://archimedes.cheme.cmu.edu/?q=nick) |
| License | Free for Problems Up to 10 Variables |
| Home Page | [The Optimization Firm](http://www.minlp.com/) |
| Download Page | [BARON Download Page](http://www.minlp.com/download) |
| MEX Interface | Available from the Download Page |
| Pre-Requisites | None |
| Version Interfaced | 12.7.7 |

### Interfacing BARON
OPTI will work with any version of BARON v12 or higher, but requires an additional interface to work with MATLAB. The Optimization Firm provides BARON free for problems up to 10 variables, while the MATLAB interface is free for any size. For larger problems contact  [Nick Sahinidis](http://www.minlp.com/contact-us) .

To get BARON, complete the following steps:

1. Visit the [Optimization Firm Download Page](http://www.minlp.com/download) and download *both* the BARON executable *and* the MATLAB/BARON interface.

1. Run the BARON installer and take note of where you installed it to. The default location is `C:\baron`.

1. From the BARON installation folder copy the BARON executable (`baron`), OpenMP library (`libiomp5md.dll`) and p-threads library (`pthreadVC2-tof.dll`) to the Interface folder within the matbar directory, i.e. matbar/Interface/ *and* rename the executable `barin`.

1. Next, unzip the MATLAB/BARON interface folder (matbar) to a suitable permanent location. 

1. Navigate to the matbar directory within MATLAB, then run `BARON_install` to install the interface.

Type the following at the MATLAB command prompt to ensure everything is working:

```matlab
>> baron(@(x) -x(1),[],[],[],0,1)
```
