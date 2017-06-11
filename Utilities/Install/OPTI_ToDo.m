%% BUG FIXES

% cbc problems lift and project cuts
% https://projects.coin-or.org/Cbc/ticket/135

% LEVMAR single inequality problem
% OPTI problems\ParameterEstimationProblem

% MOSEK SDP support


%% Ideas of functionality to add to OPTI Toolbox.
% The below will be added as time allows. Alternatively if you have an
% excellent student looking for a suitable project, I would be happy to
% co-supervise development! Obviously not all topics below are suitable,
% but definitely some for both under-grads and grads.

%% Differentiation
% - Sparse Numerical Differentiation
% Expand mklJac to provide sparse approximations of 1st derivatives,
% including making use of supplied derivative structures
%
% - Sparse + Second Derivatives Automatic Differentiation
% Expand or replace autoJac to provide exact first and second derivatives,
% including dense and sparse and automatic structure identification

%% New Problem Types
% - Bilinear Matrix Inequalities
% New constraint type

% - Minimax
% Reformulation

% - Goal Attainment
% Reformulation

% - Complementarity Problems
% New constraint type

% - Linear Least Squares
% Reformulation (QP)

% - Linear Objective + Nonlinear Constraints
% Automatic gradient

% - Quadratic Objective + Nonlinear Constraints
% Automatic hessian + gradient

% - Network Problems
% Quite a bit of work I suspect

%% Algorithm & Functionality Development
% - Automatic Piece-Wise Linear Approximations
% Automatic linearization of problems with separable nonlinear constraints
% and/or objective

% - Multi-start global optimization solver
% Use existing literature to develop a global multi-start solver which uses
% local solvers such as IPOPT within OPTI

% - Expansion of the SymBuilder framework
% Create a built in Symbolic Engine (rather than use the Symbolic Toolbox)
% plus formalize many of the methods for easier use

% - Expand into Simulink, provide optimization capabilities within Simulink
% Optimizer blocks?

% - Optimal control and dynamic optimization
% Open source dynamic optimization platform for MATLAB?

% - Formalize distribution testing
% Collect together a suite of test problems for ensuring each OPTI release
% is stable, including checking results / infeasiblities / etc

% - Improve Test Problem Suite
% The current selection of OPTI Test Problems is quite pitiful. A large
% number should be added (perhaps as an optional download?) of various
% problem types.

% - Interface to CUTEr
% For testing and checking solvers and interfaces

% - Tuning of SCIP Interface
% Hand optimize the SCIP interface for parsing nonlinear functions into
% instruction lists. Provide better support for sparse problems and other
% large-scale functions. Preallocate memory, better error checking, and so
% forth!

% - AMPL Interface update
% Preserve objective and constraint naming + other tidy ups. Perhaps
% conversion to m-file?

% - OPTI .NET Interface
% Develop a C++ .NET interface to OPTI solvers for easier solver
% deployment.

% - OPTI SCILAB Interface
% Develop a SCILAB equivalent interface for OPTI.

% - Standardize MEX Interfaces
% Many of the MEX interfaces were written ad-hoc as I had time. Many could
% do with a good clean up, and standardized via C++ classes and objects.

% - Improve plot(opti)
% Add 1D plotting capabilities, develop an algorithm to better
% automatically scale the plot, add more options for colours, lines, etc.

% - OPTITOOL
% OPTI optimization GUI, perhaps based on optimtool

