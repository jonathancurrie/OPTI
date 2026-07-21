---
title: "MOSEK"
slug: "/solvers/mosek/"
---

## MOSEK [Interfaced]
<small></small>
MOSEK is a commercial solver for convex linear and quadratic programs, including discrete variables and quadratic constraints.

|  |  |
| --- | --- |
| Category | Commercial |
| Manager | [MOSEK](http://www.mosek.com/) |
| License | Free Academic License (1 year no limits) |
| Home Page | [MOSEK](http://www.mosek.com/) |
| Download Page | [MOSEK Downloads Page](https://www.mosek.com/downloads/) |
| MEX Interface | Supplied with MOSEK |
| Pre-Requisites | None |
| Version Interfaced | 7 |

----
### Interfacing to OPTI
MOSEK is available free for 1 year for academics and includes a powerful SDP solver which is worth trying!

To get MOSEK, complete the following steps:

1. Visit the [MOSEK Academic License Request page](https://www.mosek.com/products/academic-licenses/) and fill in your details. Your license will be emailed to you.

1. Visit the [MOSEK Download page](https://www.mosek.com/downloads/) and download the executable for your system (e.g. Windows 64 bit x86).

1. Once the download has completed, start the MOSEK installer. I install using the default settings. Once installed, you must restart your computer.

1. Your license file should come through almost instantly. As per the instructions within the email, copy your license file to the folder specified. On my PC this is C:\Users\Jonathan Currie\mosek\mosek.lic

1. Next you will need to add the MOSEK MEX file to the MATLAB path. From MATLAB click File -> Set Path (or Environment -> Set Path for MATLAB > 2012a). From the path tool, click "Add Folder", navigate to the installation directory of MOSEK, and select the toolbox/r2013aom folder as below:

![mosek path](/img/opti/mosek_path.png)

If you installed in the default directory this should be similar to:

C:\Program Files\Mosek\7\toolbox\r2013aom

It is *very important you add the above folder*. If you add the wrong folder it is likely you will break the optimization toolbox and possibly OPTI.

Click OK to select the matlab folder, click "Save" to save the path changes, then type the following at the MATLAB command prompt to ensure everything is working:

```matlab
>> mosekopt
```
