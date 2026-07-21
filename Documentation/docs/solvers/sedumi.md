---
title: "SeDuMi"
slug: "/solvers/sedumi/"
---

## Self-Dual-Minimization [Interfaced]
<small></small>
SeDuMi is the defacto standard for SDP solvers, and includes functionality for solving SOCP and rotated cone problems (not yet interfaced to OPTI). For interfacing instructions see below.

|  |  |
| --- | --- |
| Category | Open Source |
| Managers | [Dr. Imre Polik](http://imre.polik.net/) and  [Prof. Johan Löfberg](http://users.isy.liu.se/en/rt/johanl/) |
| License | [GNU General Public License v3](http://www.gnu.org/licenses/gpl.html) |
| Home Page | [SeDuMi Home Page](http://sedumi.ie.lehigh.edu/) |
| Download Page | [SeDuMi Download Page](http://sedumi.ie.lehigh.edu/index.php?option=com_docman&task=cat_view&gid=69&Itemid=76) |
| MEX Interface | Supplied with SeDuMi |
| Pre-Requisites | None |
| Version Interfaced | 1.32 |

### Interfacing SeDuMi
OPTI is only interfaced to the latest release of SeDuMi (v1.31 or higher), now available on github  [here](https://github.com/sedumi/sedumi). Simply click 'zip' and download the solver, as below:

![github](/img/opti/github.png)

Once downloaded, unzip it on your computer to a suitable permanent location, and add the path to MATLAB using the following commands (make sure MATLAB's current directory is the sedumi folder):

```matlab
>> genpath(cd); 
>> savepath;
```

### Citation
J.F. Sturm, "Using SeDuMi 1.02, a MATLAB toolbox for Optimization over Symmetric Cones," *Optimization Methods and Software* 11(12), pp. 625-653, 1999
