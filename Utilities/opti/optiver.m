function current_version = optiver
%OPTI Version Information
%  Return OPTI software version. This file also contains version update
%  information.

%   Copyright (C) 2011-2015 Jonathan Currie (I2C2)

current_version = 2.16;

%History

%2.16 (..)
% - Minor bug fixes
% - Fixed bug when nlcon cannot accept x0 = zeros(n,1) during function testing (A. Khalid)

%2.15 (17/05/15)
% - Rebuilt all solvers against MKL v11.2 R3 and Visual Studio 2013
% - Added the ability to specify options in SCIP
% - Changed checkSolver to optiSolver, same API
% - Updated CLP to 1.16.6
% - Updated CBC to 2.9.4
% - Updated IPOPT to 3.12.3
% - Updated BONMIN to 1.8.2
% - Updated SCIP to 3.1.1 (note more sensitive to unbounded problems now)
% - Updated RMathlib to 3.2.0
% - Updated CppAD source to 20150429
% - Updated Documentation, mostly moved to Wiki now, tidied paths a little bit
% - OPTI source (but not binaries) is now on GitHub: https://github.com/jonathancurrie/OPTI
% - Fixed bug in convBonmin when supplying options from another solver (G.Tsoukalas)
% - Fixed bug in >=R2014b when plotting OPTI plots with bounds
% - Fixed bug in >=R2014b when plotting quadratic constraints
% - Fixed bug with bintprog no longer available in R2014b
% - Fixed bug with orange line in data fitting plots in R2014b

%2.12 (23/03/15)
% - Added warning to problems which did not solve correctly.
% - Updated CLP to v1.16.5, removed Aboca support
% - Updated CBC to v2.9.3
% - Updated IPOPT to v3.12.1
% - Updated BONMIN to 1.8.1
% - Re-added CLP to be able to solve QPs (earlier bug appears fixed)
% - Rebuilt all solvers against MKL v11.2 R2

%2.11 (16/09/14)
% - Rebuilt all solvers against MKL v11.2 R0
% - Added SymBuilder Examples page
% - Added extra error information to SCIP function processing
% - Fixed bug in SCNLE sparse test (Stephen)
% - Fixed bug in confidence bound plotting with unsorted xdata
% - Fixed bug in DNLS problems with initialT ignored for single sampling rate problems

%2.10 (20/07/14)
% - Reorganised OPTI folder structure to reduce number of folders
% - Removed support for lipsol and qsopt
% - Removed demos, re-written as examples
% - Rebuilt all solvers against MKL v11.1 R3
% - Added support for intlinprog (2014a) (including opti_ overload)
% - Added ability to calculate fit statistics for NLS/DNLS problems
% - Added fitting function 'optifit'
% - Added experimental symbolic derivative & solving function 'optisym'
% - Added automatic plotting of prediction confidence bounds
% - Added RMathlib with a library of statistical functions
% - Added code generation (including automatic derivatives via CppAD) from SymBuilder models
% - Enabled MOSEK install test
% - Fixed options bug in opti_nomad with Inf bounds (M. Porcelli)
% - Fixed bug in nlopt displaying suboptimizer when not used
% - Fixed bug in SymBuilder with 0.5 in quadratic programs
% - Fixed bug in Cbc incorrectly missing infeasibility in some models
% - Fixed bug in derivative checker with Hessian with no NL constraints
% - Updated IPOPT to v3.11.8
% - Updated NLOPT to v2.4.2
% - Data fitting weights are now treated as sqrt(weights) internally 

%2.05 (29/12/13)
% - Rebuilt all solvers against MKL v11.1 R1
% - Added support for BARON (Global MINLP solver)
% - (Re)Added support for MOSEK ((MI)LP/QP/QCQP solver)
% - Added ability to solve SNLE/SCNLE problems as NLP without dummy objective
% - Added ability to plot SLE/SNLE/SCNLE problems
% - Added ability to supply weights for NLS and DNLS problems
% - Added ability to change fmincon algorithm via optimset options
% - Added plotting of the initial guess to plot(optiObj)
% - Added build date to MEX files
% - Added ability to specify 'solver' option in opti arg list
% - Added ability to read mathematical files directly into opti constructor
% - Added cstepJac, cstepHess and cstepHessLag for complex step first and second derivatives
% - Added mkltrnlsset for setting MKLTRNLS options
% - Updated Clp to v1.15.6 (removed QP support due to heap corruption)
% - Updated Cbc to v2.8.8
% - Updated Bonmin to v1.7.4
% - Updated SCIP to v3.0.2
% - Updated IPOPT to v3.11.7, added PARDISO as a linear solver option
% - Updated NOMAD to v3.6.2
% - Updated L-BFGS-B interface, added lbfgsbset
% - Updated CBC interface, added presolver + cbcset, changed return codes + arg order
% - Updated support for SeDuMi v1.32 and higher
% - Smaller fonts on higher dimensional plots
% - Fixed slow code for checking integer variable arrays
% - Fixed bug in Clp when suppling H with no nonzero terms (J.Lofberg)
% - Fixed bug when writing MPS files with SOS
% - Fixed bug in HYBRJ, LMDER and LEVMAR with large problems
% - Fixed bug in symJac, symHess, symHessLag with scalar problems
% - Fixed bug in mklJac when supplied sparse data (now returns error)

%v2.00 (20/06/13)
% - Rebuilt all solvers against R2013a, VC++ 2012 and MKL v11 R4
% - Updated pre-requisite checking for VC++ 2012 and Intel C++ XE 2013
% - Updated SeDuMi interface to support v1.31 and above only
% - Updated amplRead to return quadratic row constraints
% - Updated NOMAD to accept problems up to 5000 variables
% - Updated 'set' routines to create common code checks
% - Updated OOQP MEX interface, added MA57 as a linear solver option (also
%   now tril H and no A' as well as Ctrl C, maxtime and maxiter options) 
% - Updated symJac, added symHess, symHessLag, also returns
%   sparsity pattern + more robust and added symDynJac for DNLS problems
% - Updated BONMIN to v1.7.0
% - Updated Clp to v1.15.1 (+ changes to argument order + added presolver + clpset + Aboca)
% - Updated CBC to v2.8.2
% - Updated CoinUtils to v2.9.1
% - Updated GLPK to v4.48
% - Updated NOMAD to v3.6.1 (+ changes to nomadset)
% - Updated IPOPT to v3.11.1
% - Updated CppAD in SCIP to v20130422
% - Updated ASL to v20130419
% - Added LIPSOL (LP solver) 
% - Added basic dynamic parameter estimation solving (+ optidynset)
% - Added basic multi-start solver (multisolve)
% - Added automatic derivative checker (optiset option 'derivCheck')
% - Added most remaining IPOPT options to ipoptset
% - Added ability to write GAMS files from OPTI (inc NL via SCIP)
% - Added automatic VS project builder code to most solvers
% - Added objbias term for linear and quadratic problems (inc File IO)
% - Added x0 as an input argument to opti_linprog and opti_quadprog
% - Added lambda as an output argument from opti_linprog and opti_quadprog
% - Added complex MUMPS build
% - Added ability for IPOPT, NOMAD, NLOPT to solve SNLE, NLS
% - Added ability to plot bounded or x0 supplied unsolved problems
% - Added ability to plot 1D to 5D problems
% - Fixed bug in NOMAD with >= nonlinear constraints
% - Fixed bug in IPOPT when not supplying opts.ipopt
% - Fixed bug in BONMIN when not supplying any constraints
% - Fixed bug in BONMIN when supplying double sided row constraints
% - Fixed bug in sdpWrite with row constraints
% - Fixed bug in SCIP with quadratic constraints (incorrect l index)
% - Fixed bug in optiprob with ndec that overwrites x0
% - Removed support for CPLEX v12.5.0.0 with R2013a x64 (strange crashes)

%v1.81 (02/04/13)
% - Added ability to specify quadratic constraints in row form
% - Added support for plotting quadratic equalities
% - Changed how SOS are stored internally
% - Changed xtype specification for CBC
% - Changed SOS specification for CBC, LP_SOLVE and SCIP and coinW
% - Rebuilt CSDP with aligned memory allocation
% - Rebuilt solvers and utilities against MKL 11.0 R3
% - Fixed bug in filterSD glcpd.f
% - Fixed bug in filterSD [sparse] with workspace memory allocation
% - Fixed bug in filterSD [sparse] with bounded problems
% - Fixed bug in structure generation for linear constraints with nonlinear solvers
% - Fixed bug in MIQCQPs with CPLEX

%v1.80 (27/03/13) 
% - Added Semidefinite Programming Constraints (sdcone, sedumi)
% - Added Semidefinite cone plotting
% - Added DSDP (SDP solver)
% - Added CSDP (SDP solver, parallelized with OpenMP)
% - Added interface to SeDuMi
% - Added ability to read SeDuMi arguments into opti()
% - Added sdpRead() to read SDPA & SeDuMi problems into opti()
% - Added sdpWrite() to write OPTI problems to SDPA & SeDuMi files
% - Added some of the SDPLIB problems to opti test problems
% - Changed warning levels to 'all', 'critical' (default) or 'none'
% - Changed checkSolver to also provide solver configuration details
% - Changed IPOPT & BONMIN to accept linear constraints directly

%v1.79 (05/03/13)
% - Substantial updates to ASL mex interface, now reads QCQPs natively
% - Updated amplRead to allow linear constraints to be identified for NLPs
% - AMPL constraint linearity now available via field 'conlin' returned from amplRead
% - Changed solve(opti) to automatically write AMPL .sol files once an AMPL problem is solved
% - Changed CPLEX OPTI interface to use CPLEX class object
% - Rewrote and tidied up NLOPT MEX interface
% - Added #FuncEvals returned from NLOPT
% - Fixed bug in buildOpti QCQP checking of l
% - Fixed bug in asl when requesting Hessian structure multiple times
% - Fixed bug when plotting semidefinite / indefinite QCs
% - Fixed bug in SCIP QCQPs
% - Fixed bug in NLOPT Jacobian callback

%v1.78 (17/02/13)
% - Updated SCIP to correct for SOC bug (also updated CppAD to 20130205)
% - Updated SCIP to allow reading and solving of AMPL .nl models.

%v1.77 (04/02/13)
% - Changed error checking to enforce double precision numerical arguments only
% - Fixed IPOPT iteration callback to return complete x vector
% - Rebuilt solvers and utilities against R2012b and MKL 11.0 R2
% - Added To-Do list for possible contributors to OPTI (see Help folder)

%v1.76 (26/01/13)
% - Added filterSD (NLP solver)
% - Updated amplRead to interface to ampl.exe directly
% - Updated SCIP to v3.0.1 and SoPlex to v1.7.1 (unfortunately doesn't fix 32bit crashes)
% - Minor documentation updates
% - Minor bug fixes

%v1.75 (20/01/13)
% - Updated SCIP interface to solve Global NLPs and MINLPs (beta)
% - Updated BONMIN interface to allow use of CPLEX as an MILP solver
% - Changed error checking to enforce dense problem vectors
% - Updates to checkSol()
% - Changes to allow IPOPT and BONMIN to solve (MI)QCQPs via solver choice
% - Added extra return argument structure to IPOPT iteration callback
% - Minor updates to interfaces for LBFGSB, MKLTRNLS, HYBRJ and NL2SOL solver
% - Updates to ipoptset and bonminset

%v1.72 (24/12/12)
% - Changed IPOPT return function evaluation names to match ipopt.m description
% - Minor documentation updates
% - Updated CLP to v1.14.8
% - Updated CoinUtils to v2.8.8
% - Updated CBC to v2.7.8

%v1.71 (13/11/12)
% - Added ability for IPOPT (and BONMIN & SCIP) to use Matlab's supplied MA57 linear solver
% - Rebuilt solvers and utilities against Intel MKL 11.0 R1
% - Fixed bug in opti_fsolve for scalar problems

%v1.70 (14/09/12)
% - Added SCIP v3.0.0
% - Added M1QN3 v3.3
% - Fixed memory leak in OOQP interface
% - Fixed matrix ordering problem in coinRead
% - Rebuilt solvers and utilities against Intel MKL 11.0 R0
% - Updated CLP to v1.14.7
% - Updated CoinUtils to v2.8.7
% - Updated LEVMAR to v2.6
% - Updated NLOPT to v2.3
% - Updated NL2SOL to v2.3, now handles bounds
% - Updated BONMIN to latest stable build
% - Migrated solver documentation to Wiki
% - Updates to installer

%v1.64 (25/08/12)
% - Fixed documentation bug with sense (1 is min, -1 is max)
% - Updated ASL interface to v20120624
% - Fixed ASL interface bug with undefined operations

%v1.63 (07/08/12)
% - Fixed infeasible problem status bugs in nomad & ooqp
% - Removed support for MOSEK

%v1.62 (11/07/12)
% - Rebuilt solvers and utilities against Intel MKL 10.3 R11
% - Updated solvers to return #FuncEvals where possible
% - Changed info field 'StatusString' to 'Status'
% - Updates to SymBuilder

%v1.60 (18/06/12)
% - Updated error catching on MATLAB callbacks
% - Updated BONMIN to v1.6.0
% - Updated CBC to v2.7.7
% - Fixed bug in autoJac with atan
% - Added couple of small utilities (SymBuilder + VisualStudio Builder)

%v1.58 (20/05/12)
% - Tidied up Hessian specification (few bugs here)
% - Updated fmincon interface to use Hessian
% - Updated amplRead to add Hessian information by default
% - Updated demos to reflect Hessian changes
% - Rebuilt all MEX files against R2012a
% - Removed Hmul from problem description (unlikely to be used)

%v1.57 (14/05/12)
% - Fixed memory bug in amplRead QP reader

%v1.56 (11/05/12)
% - Added optiset option iterfun for iteration callbacks
% - Updated all NLP and NLS solvers to utilize iteration callbacks
% - Removed solverset functions including solver name

%v1.55 (04/05/12)
% - Added optiset option maxfeval
% - Added functionality to solve MIQP / QCQP / MIQCQP as NLP / MINLP
% - Changed opti constructor to now allow problem parameters directly (remains backwards compatible)
% - Changed optiset option tolfun to tolrfun and tolafun
% - Changed solverset functions to include solver name for option checking
% - Changed optiprob display format + added probtype to override OPTI ptype decision
% - Updated all examples and documentation to the new opti constructor format
% - Updated nlopt interface to include display + ctrl c exit
% - Updated all compatible solvers with new maxfeval and tol options
% - Updated checkSolver to include a 'config' option for checking options
% - Updated many nonlinear solver opti_* interfaces for easier user access
% - Updated NLS plotting routine
% - Updated warning messages to avoid stacktrace
% - Removed OPTI mfile solvers (not competitive)
% - Removed convert() method

%v1.50 (26/04/12)
% - Added NOMAD Global NLP solver
% - Added ability to specify linear constraints in row form (rl <= Ax <= ru)
% - Added ability to specify nonlinear constraints in row form (cl <= c(x) <= cu)
% - Added Global NLP demo
% - Updated clp, cbc, ooqp, ipopt to use row linear constraints
% - Updated coinRead, coinWrite, amplRead to use row linear constraints
% - Updated ipopt to use row nonlinear constraints
% - Updated amplRead to use row nonlinear constraints
% - Updated most constrained solvers to return lambda (Lagrange multipliers at solution)
% - Updated pswarm to enable vectorized function calls
% - Added interface to MATLAB Global Optimization Toolbox
% - Fixed constraint gradient bug in MATLAB fmincon interface
% - Changed checkSolver to require 'nls' instead of 'all_nls'
% - Updates to uninstaller
% - Small bug fixes

%v1.34 (14/03/12)
% - Changed install test to check relative error
% - Updated checkSolver to display solver version and constraint information

%v1.33 (13/03/12)
% - Added PSwarm Global Optimization Solver
% - Added interface to MOSEK
% - Added interfaced solver installation guide
% - Fixed bug with row LHS linear equality constraints

%v1.31 (02/03/12)
% - Added Intel MKL Trust Region Bounded NLS Solver

%v1.30 (01/03/12)
% - Updated BONMIN to v1.5.2
% - Updated IPOPT to v3.10.2
% - Updated GLPK to v4.47
% - Updated CBC to v2.7.6
% - Updated CLP to v1.14.6
% - Updated CoinUtils to v2.8.6 (File IO)
% - Fixed bug with MINLP problems with no bounds specified
% - Fixed bug with GMPL reader
% - Added objective and gradient pre-checks to IPOPT interface

%v1.28 (22/02/12)
% - Fixed bug with checkSolver locating .m instead of .mex

%v1.27 (03/02/12)
% - Fixed bug with HELP location

%v1.26 (18/12/11)
% - Fixed bug with SOS constraints in CBC + LPSOLVE wrappers

%v1.25 (15/12/11)
% - Minor updates to MINPACK + NL2SOL solvers
% - Added Matlab Overloads demo
% - Added Basic Functionality demo

% v1.24 (12/12/11)
% - Setup NL2SOL + LM_DER to solve SNLE problems
% - Added OS check
% - Added VC++ 2010 pre-req check

% v1.23 (09/12/11)
% - Added SNLE problems
% - Added HYBRJ
% - Added Special Ordered Sets (SOS)
% - Added SOS handling to CPLEX, CBC & LP_SOLVE interfaces
% - Added opti_mintprog
% - Upgraded File IO (MPS) to read / write SOS
% - Upgraded CLP to v1.14.5
% - Upgraded CBC to v2.7.5
% - Upgraded CoinUtils to v2.8.5 (File IO)

% v1.20 (11/11/11)
% - Added Matlab Optimization Toolbox Overloads
% - Added m-file help documentation to MEX functions
% - MEX files rebuilt in 2011b

% v1.15b (19/10/11)
% - OPTI Toolbox released under the BSD 3 Clause License

% v1.15 (15/10/11)
% - Added NL2SOL
% - Added amplRead
% - Added AMPL documentation
% - Updated display
% - Fixed opti_miqp fval bug

% v1.13 (05/10/11)
% - Added MINPACK LM_DER

% v1.12 (03/10/11)
% - Added NLS problems
% - Added LEVMAR
% - Updated display

% v1.10 (12/09/11)
% - Fixed bug in opti_minlp with split integer variables
% - Added BONMIN
% - Added OOQP built against Intel MKL PARDISO
% - Added L-BFGS-B

% v1.05 (31/08/11)
% - Added adiff for automatic differentiation, re-organised diff folder
% - Added differentiation demo
% - Fixed bug in uninstaller removing control toolbox options path (oops)

% v1.00 (25/08/11)
% - Added hatching to quadratic and nonlinear constraints
% - Fixed display iterations requiring matlab drawnow
% - Added Ctrl-C detection to compatible solvers
% - Added CBC
% - Added Post Install Test
% - Fixed qsoptmex bug editing Matlab bounds

% v0.95 (20/08/11)
% - Identified bug in R2010a/b with opti_mumps fortran error. Tried
% rebuilding from R2010a still crashes, tried rebuilding using Ifort 2011
% release 5, still crashes, not sure here...
% - Added MKL djacobi for nonlinear gradient and Jacobian estimation
% - Completely overhauled options setup, now more flexible for all solvers, 
% and more powerful for NLP solvers.
% - Added convert() to convert between NLP problem formats
% - Added coinRead to read MPS, QPS, GMPL, GAMS and LP problem formats
% using CoinUtils.
% - Added coinWrite to write MPS, QPS and LP problem formats using
% CoinUtils.
% - Formalized SLE problems
% - Added MIQCQP problems
% - Templated SNLE + SCNLE problems
% - Added NL plotting
% - Added UNO problems
% - Updated HTML documentation
% - Overhauled plotting of problems, now a common file controls overlays
% - Fixed bug(s) in linear constraint drawing
% - Changed all examples to demos, added to documentation
% - Added linear -> nonlinear constraint conversion for NLPs
% - Added automatic infeasible x0 moving for bounded NLPs with NLOPT
% - Finally finished opti_minlp!

% v0.9 (05/08/11)
% - Another pre-release
% - Updates and bug fixes to opti_readMPS
% - Added MIQP problems + OPTI MIQP solver
% - Added plotting of MIQP problems
% - Updated clpmex.cpp to handle sparse A + Quadratic Objective
% - Added QSopt
% - Added html documentation (under development)
% - Re-organised utilities folder structure
% - Updated CLP to v1.14.0 + build info
% - Updated IPOPT to v3.10.0 + build info

% v0.85 (27/07/11)
% - Another pre-release
% - Added nlp_HS to auto generate HS NLP problems
% - Added test_Solver for individual solver benchmarking
% - Added getSymJac for Symbolic generation of Jacobian
% - Added LP,QP,MILP,NLP examples
% - Re-wrote optiBench
% - Added extra optiset options (under development)
% - Added opti_readMPS to read MPS/QPS files (under development)

% v0.81 (24/07/11)
% - Another pre-release
% - Small updates to display
% - Removed need for numdif if derivative free optimizer is used

% v0.8 (22/07/11)
% - Another pre-release
% - Added NLOPT (still under development)
% - Removed all .lib files from distribution (licensing problems + size)
% - Fixed opti_ipopt bug with constraint bounds
% - Added display function to opti object (under development)

% v0.7 (17/07/11)
% - Another pre-release
% - Added IPOPT + OOQP
% - Rewrote checkSolver
% - Updated optiBench
% - Added convMatlab to return fmincon structure
% - Added NLP solving functionality
% - Rewrote optiset

% v0.5 (11/07/11)
% - Another pre-release
% - Added CLP + MUMPS
% - Added detailed solver compiling instructions
% - OPTI milp solver re-written

% v0.3 (08/07/11)
% - Another pre-release
% - Added GLPK + LP_Solve
% - Added optiobj and opticon
% - Re-organised folder structure

% v0.1 (23/06/11)
% - Pre-release of OPTI Toolbox for DIW.

end

