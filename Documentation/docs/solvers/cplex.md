---
title: "CPLEX"
slug: "/solvers/cplex/"
---

## CPLEX [Interfaced]
<small></small>
CPLEX is a commercial solver for linear and quadratic programs, including discrete variables and quadratic constraints.

|  |  |
| --- | --- |
| Category | Commercial |
| Manager | [IBM ILOG](http://www-01.ibm.com/software/websphere/products/optimization/) |
| License | Free Academic License |
| Home Page | [CPLEX](http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/) |
| Download Page | [IBM Academic Initiative](https://www.ibm.com/support/pages/node/138207) |
| MEX Interface | Supplied with CPLEX |
| Pre-Requisites | None |
| Version Interfaced | 12.6 |

### No 64bit Support
Since Cplex v12.5 there is a bug in the Cplex class object under MATLAB 64bit which can crash MATLAB. I raised this issue with IBM/ILOG in July 2013 however I have yet to obtain a solution. Therefore 64bit support of Cplex under OPTI is currently disabled.

### Interfacing CPLEX
OPTI will work with any version of Cplex v12.3 or higher. IBM now very kindly issues CPLEX license free for academics, and is a fantastic solver so it is well worth interfacing!

To get CPLEX, complete the following steps:

1. Visit the [IBM ILOG Academic Initiative page](http://www-03.ibm.com/ibm/university/academic/pub/page/mem_join) and register for the IBM Academic Initiative.

1. Visit [IBM Software Catalog page](https://www.ibm.com/support/pages/node/138207), where you will need to log in to your IBM account.

1. Enter "cplex" in the "Find by search text" textbox. Expand the node next to the latest version of IBM ILOG CPLEX Optimization Studio and select the tickbox next to "... for Windows x86-64 Multilingual". Alternatively select the x86-32 box if you want the 32 bit version.

1. Read and agree to the license conditions of the download, click I Agree, and then click Download now.

1. Once the download has completed, start the CPLEX installer. You may need to run in compatibility mode for Windows 7 if you are running Windows 8. I install to the default location.

1. Once installed, you will need to add the CPLEX Toolbox to the MATLAB path. From MATLAB click File -> Set Path (or Environment -> Set Path for MATLAB > 2012a). From the path tool, click "Add with Subfolders", navigate to the installation directory of CPLEX, and select the toolbox folder as below:

![cplex path](/img/opti/cplex_path.png)

If you installed in the default directory this should be similar to:

C:\Program Files\IBM\ILOG\CPLEX_Studio126\cplex\matlab

Click OK to select the matlab folder, click "Save" to save the path changes, then type the following at the MATLAB command prompt to ensure everything is working:

```matlab
>> Cplex
```
