Copyright (C) 2009 A. Ismael F. Vaz and L. N. Vicente.

http://www.norg.uminho.pt/aivaz http://www.mat.uc.pt/~lnv

http://www.norg.uminho.pt/aivaz/pswarm

This library is free software; you can redistribute it and/or modify
it  under  the  terms  of  the  GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,  but
WITHOUT  ANY  WARRANTY;  without  even  the  implied   warranty   of
MERCHANTABILITY  or  FITNESS  FOR  A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the  GNU  Lesser  General  Public
License  along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
USA.



                     ===========================
                     Parallel PSwarm Version 1.5
                     ===========================

                           ---------------
                             Read me file
                           ---------------


Welcome to Parallel Pattern Swarm version 1.5 (PPSwarm v1.5)

This code provides an implementation of the algorithm described in:

A.I.F.  Vaz  and L.N.Vicente, A particle swarm pattern search method
for   bound  constrained  global  optimization,  Journal  of  Global
Optimization, 39 (2007) 197-219.

The linear constrained approach is developed in the paper:

A.I.F.  Vaz  and  L.N. Vicente. PSwarm: A hybrid solver for linearly
constrained global derivative-free optimization.


Available online at:

www.norg.uminho.pt/aivaz/pswarm


This file describes the following topics:


1. Change log
2. System requirements
3. Compilation/Installation
4. Contents of the ppswarm_v1_1 directory
5. The command line
6. Credits


1. Change log
----------------------------------------
PSwarm v0.1 implements the algorithm described in

A.  I.  F.  Vaz  and  L. N. Vicente, A particle swarm pattern search
method  for bound constrained global optimization, Journal of Global
Optimization,    39    (2007)   197-219.   (available   online   at:
http://dx.doi.org/10.1007/s10898-007-9133-5)

PSwarm  v0.1  provides  a  standalone  version  (where the user must
provide  a simple bound optimization problem) that can run in serial
or parallel mode.

PSwarm  v0.1  provides also an AMPL (www.ampl.com) interface version
in serial mode.

Version v1.1 implements the algorithm described in

A.I.F.  Vaz  and  L.N. Vicente. PSwarm: A hybrid solver for linearly
constrained global derivative-free optimization.

PSwarm   v1.1   is  now  able  to  handles  both  bound  and  linear
constraints.

PSwarm   v1.1   provides   a   R   (www.r-project.org)   and  Python
(www.python.org) interface in addition to the AMPL interface.

PSwarm  v1.2  corrects  a  bug  in  the  pswarm_py.c file and in the
makefile.

PSwarm v1.3 allows user to provide a iteration print function called
each  iter  mod  iprint  = 0, where iter is the iteration number and
iprint  is  a  new  option available in the algorithm. These options
(with many others) are now exported to R and Python.

PSwarm v1.4 calls the objective function in a vectorized way.

PSwarm v1.5   corrects   some  memory  allocation,  allowing  bigger
problems to be addressed.


2. System requirements
----------------------------------------

PSwarm  v1.4  can  be  compiled  with  or without linear constraints
support.

This option is available at compilation time allowing users to still
use the PSwarm algorithm in different systems configurations.


PSwarm with only bound constraints support:
 -  C  compiler  (gcc  for  the  serial Linux version, mpicc for the
 parallel Linux version or VC++ for the serial Windows version);
 - AMPL interface library for the AMPL version.

PSwarm with bound and linear constraints:
 -  C  compiler  (gcc  for  the  serial Linux version, mpicc for the
 parallel Linux version or VC++ for the serial Windows version);
 -  BLAS  (Basic Linear Algebra Subroutines) library (may be already
   available  on your system, otherwise download it from the web and
   make a local BLAS library);
 - LAPACK (Linear Algebra PaCKage) library (may be already available
   on  your  system,  otherwise  download it from the web and make a
   local LAPACK library);
 - AMPL interface library for the AMPL version;
 - A system Python installation for the Python dll version;
 - A system R installation for the R dll version.


3. Installation
----------------------------------------

The code requires  approximately  12.5Mb  of  hard  disk  space
(source files, windows libraries, AMPL student edition binaries,
etc).

Unzip  the PPSwarm_v1_1.zip file in any location using your favorite
unzip  software. A directory PPSwarm_v1_1 will be created, having as
subdirectories  'nl',  'include',  'libs',  which refers to the AMPL
problems  models  (in  the  .nl format), the corresponding interface
include  files  and a directory to place additional libs (AMPL, BLAS
or LAPACK), respectively.

Two makefiles are provided to  assist  in  the  solver  compilation.
'makefile.vc' is for the Microsoft Windows Visual C++ and 'makefile'
is for Linux.


LINUX:
------
A makefile for Linux is provided.

In  order to obtain the PSwarm with the linear constraints support a
'_linear'  suffix  should  be  added  to the make command line (e.g.
'make  ampl_linear' will make the pswarm binary with AMPL and linear
constraints support).

Type  'make  serial'  for  the  standalone  serial  version.  In the
standalone  version  the user must provide the problem definition in
the  file  user.c. The 'objfun' C function is expected to return the
objective  function value at 'x'. The 'set_problem_dimension' return
the number of problem variables and 'set_problem' sets the lower and
upper  bounds  on  the variables together with an initial guess. The
user may use the 'user_init' function to initialize some data.

Type 'make parallel' for a standalone parallel version (see previous
comments  on  how  to tune the objective function and bounds). There
is  no manual for the parallel version, but the code just assigns to
each  processor  the  computation of an objective function value. We
advice   to  have  2*(n+1)+1  processor  (since  currently  we  have
2*(n+1)   directions  in  the  poll  step  plus  one  processor  for
running  the  algorithm). See file 'job.pbs' for an example on using
pswarm in a PBS system.

Type 'make ampl' to make the ampl version of pswarm.

Type 'make lib' to make the a pswarm library (libpswarm.a). The user
functions defined in user.c are not included and must be provided by
the user.

Type  'make  r' to make the pswarm dynamic library to be used with R
(see the provided .r example).

Type  'make  py'  to make the pswarm dynamic library to be used with
Python (see the provided .py example).

WINDOWS:
--------
A makefile.vc for Microsoft VC++ is provided. The same comment
applies to the _linear suffix in the make command line arguments.

Type 'nmake -f makefile.vc serial' for the serial standalone
version.

Type 'nmake -f makefile.vc ampl' for the AMPL version.

Type 'nmake -f makefile.vc rr' for the R dynamic library. Note that
under Windows its 'rr' for the version without linear constraints
and 'r_linear' for the version with linear constraints.

Type 'nmake -f makefile.vc py' for the Python dynamic library.

Type 'nmake -f makefile.vc lib' for the pswarm library.


WINDOWS MPI version is not provided.


4. Contents of the PPSwarm_v1_5 directory
----------------------------------------

In the directory PPSwarm_v1_5 there are the following files:

------ README -------------
README.txt - This readme file

lgpl.txt - GNU LESSER GENERAL PUBLIC LICENSE


------ Solver files -------
pswarm_main.c(h) - Main file for the solver (interface with AMPL and
                   problem definition).

pswarm.c(h) - Main algorithm (the particle swarm).

pattern.c(h) - Pattern search algorithm (the poll step).

user.c - User provided functions for the standalone version.

cache.c(h) - implements a cache for the objective function. Disabled
             by  default.  If  this  feature  is  to be enabled then
             comment  out  the  load_cache_file  and save_cache_file
             routines  call  in  pswarm_main.c.  In  your  objective
             function  a  call  to the cache should also be provided.
             This is an untested code.

mve_presolve.c  -  Interior  point  code to compute a point interior
                 to  the  ellipsoid  (initial guess to the ellipsoid
                 center).

mve_solve.c  -  Computes the Maximum Volume Ellipsoid by an interior
              point code.

pswarm_py.c(h) - Python interface to PSwarm.

pswarm_r.c(h) - R interface to PSwarm.

RunPswarm.py  -  A  Python script to illustrate the PSwarm use under
               Python.

RunPswarm.r - A R script to illustrate the PSwarm use under R.

hs024.py  -  The Hock and Schittkowski problem description in Python
           (used by RunPSwarm.py)

hs024.r  -  The Hock and Schittkowski problem description in R (used
          by RunPSwarm.r).


------ Directories --------
nl - AMPL models (in .nl format).

include - AMPL (or additional) include files needed for the
          interface.

libs   -  Directory  to  place  additional  libraries  (AMPL,  BLAS,
          LAPACK).

------ Miscellaneous ------
job.pbs - a script for the PBS system.

showcache.c - display the cache file contents.

makefile - The Linux makefile.

makefile.vc - The Windows makefile.



5. The command line
----------------------------------------

For the standalone  version  just  type  'pswarm'  to  run  for  the
compiled problem.

For  the  AMPL  version  use 'pswarm stub.nl' or use the ampl binary
to solve the problems provided in the AMPL language.

For the parallel version see the 'job.pbs' file for a PBS script.

To run the Python example just type 'from RunPSwarm import *' in the
python command line.

To  run  the  R  example  just type 'source('RunPSwarm.r')' in the R
command line.


6. Credits
----------------------------------------

The MVE code is a translation into C of the MVE MATLAB interior
point code developed by Zhang and Gao:

Y.~Zhang  and  L.~Gao,  "On numerical solution of the maximum volume
ellipsoid  problem",  SIAM  Journal  on  Optimization, 14(1):53--76,
2003. http://www.caam.rice.edu/~zhang/mve/index.html

The  Windows  BLAS and LAPACK libraries can be obtained from Karl M.
Syring at http://www.weihenstephan.de/~syring/f2c/f2c.html
