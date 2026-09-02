---
title: "FAQ"
slug: "/project/faq/"
---

Below are some of the common questions I get from OPTI users. Please browse here before emailing me for technical support!

## General OPTI Questions
### Does OPTI work on Mac or Linux?
While MATLAB code will run on any OS, the solvers supplied with OPTI have only been compiled against Windows x64 (64bit). However as all MEX interface source is supplied, as well as detailed compilation instructions, you can always have a go yourself at compiling the solvers you want on your OS! **Note I am unable to assist with compiling for Linux or Mac.**

### What pre-requisites does OPTI require?
OPTI requires MATLAB 2011a or above to work, as well as the  [Visual C++ 2017](https://go.microsoft.com/fwlink/?LinkId=746572) and  [Intel Fortran XE 2019](https://software.intel.com/en-us/articles/redistributable-libraries-for-intel-c-and-fortran-2019-compilers-for-windows) Runtimes. When installing the runtimes, ensure you choose the version that matches your MATLAB version (i.e. 64bit), and always download the latest Fortran update from Intel.

### How do I install OPTI?
The best way to get OPTI is to clone it using Git from  [GitHub](https://github.com/jonathancurrie/OPTI) . If you are unsure how to use Git, see the read me on the GitHub page. Once you have cloned/downloaded OPTI, in MATLAB 'cd' to the download directory, and run opti_Install.m. If you have a previous version of OPTI, it will be automatically removed from the MATLAB path (you will still need to delete it manually afterwards), and the new version will be ready to go! If you do not save the path changes you will need to run the OPTI installer each time you start MATLAB.

### How do I customize solver-specific options? {#s-opts}
While [optiset](../guides/advanced/opts.md) provides a few of the common optimization settings, once you get serious you will want to start customizing the solver-specific options to obtain a better/faster solution. To do this follow the example [here](../guides/advanced/opts.md#s-opts).

### Solving my problem using OPTI causes MATLAB to crash! 
Every so often you may pose a problem that either my interface (likely) or the solver (unlikely) is not set up to solve correctly. As all solvers are compiled MEX files, some exceptions and memory errors can crash MATLAB. These problems need to be reported so I can fix them, so please post a question on the  [OPTI Forum](https://groups.google.com/forum/#!forum/opti-toolbox-forum) with details of your problem, and the error. 

### I have a multi-core PC but OPTI only uses 1 core?
Unfortunately most solvers supplied with OPTI are not multi-threaded (parallelized). This means no matter how many cores you have, the solvers themselves will only use 1 core. I do try and leverage parallelism by linking against Intel MKL which provides multi-threaded BLAS and LAPACK libraries, but these are low-level and may not improve performance much. There are some solvers which are parallelized such as [CSDP](../solvers/csdp.md), and from OPTI v2.01 [IPOPT](../solvers/ipopt.md) can leverage a parallel linear solver, but only very large problems will actually benefit. Note by default, parallelization settings are disabled so check the relevant options method for each solver.

### I get the error "Maximum variable size allowed by the program is exceeded"
OPTI can solve really big problems by exploiting sparsity, with no size limits imposed on any problem within the *software*. However ultimately the maximum size is limited by your *hardware*, specifically RAM. If you are getting this error MATLAB has run out of memory to store your problem. Ensure you are using 64bit MATLAB, then it is time to buy more RAM. Remember RAM is really cheap and the easiest upgrade you can do to any PC, so it is worth considering.

### Can you add my solver to OPTI?
Sorry - OPTI is no longer in development.

### How do I cite OPTI?
Have a look at the citation page [here](./citation.md).

## Solver Operation Questions
### Can I force a solver to stop?
Most solvers are setup to use the standard Ctrl C (press and hold Ctrl and C) exit. If you are lucky, the solver may even return the last solution found! Some solvers may take a few seconds to register the key press, or you may need to do it a few times. This is due to the MEX interface manually checking for this exit only when *allowed* to by the solver.

### I've pressed Ctrl C a number of times, but it won't stop!
If your solver supports Ctrl C termination (check `optiSolver('config')`) and it still won't exit, it may be *stuck* in the MATLAB callback. Ctrl C will only exit the MEX function, it will not exit MATLAB code called from the MEX function. You can either wait - or kill the MATLAB process.

### I've set the option `maxtime` but the solver keeps going!
Check the solver supports this option (check `optiSolver('config')`).
