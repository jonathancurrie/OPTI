% SCIP  Solve a LP/MILP/QP/MIQP/QCQP/MIQCQP/NLP/MINLP using SCIP
%
% SCIP uses the Solving Constraint Integer Programs library, available
% from: http://scip.zib.de/scip.shtml
%
%   [x,fval,exitflag,stats] = scip(H, f, A, rl, ru, lb, ub, xtype, sos, qc, nl, opts)
%
%   Input arguments*:
%       H - quadratic objective matrix (sparse, optional [NOT TRIL / TRIU])
%       f - linear objective vector
%       A - linear constraint matrix (sparse)
%       rl - linear constraint lhs
%       ru - linear constraint rhs
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       xtype - string of variable integrality ('c' continuous, 'i' integer, 'b' binary)
%       sos - SOS structure with fields type, index and weight (see below)
%       qc - Quadratic Constraints structure with fields Q, l, qrl and qru (see below)
%       nl - Nonlinear Objective and Constraints structure (see below)
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       stats - statistics structure
%
%   Option Fields (all optional - also see scipset):
%       tolrfun - LP primal convergence tolerance
%       maxiter - maximum LP solver iterations
%       maxnodes - maximum nodes to explore
%       maxtime - maximum execution time [s]
%       display - solver display level [0-5]
%       objbias - constant objective bias term
%
%   Return Status:
%       0 - Unknown
%       1 - User Interrupted
%       2 - Node Limit Reached
%       3 - Total Node Limit Reached
%       4 - Stall Node Limit Reached
%       5 - Time Limit Reached
%       6 - Memory Limit Reached
%       7 - Gap Limit Reached
%       8 - Solution Limit Reached
%       9 - Solution Improvement Limit Reached
%      10 - Restart Limit Reached
%      11 - Problem Solved to Optimality
%      12 - Problem is Infeasible
%      13 - Problem is Unbounded
%      14 - Problem is Either Infeasible or Unbounded
%
%   Special Ordered Sets (SOS):
%       type    - A string containing '1' or '2' for each SOS1 or SOS2 constraint.
%       index   - A double array of the indices of the variables in the SOS,
%                 group multiple SOS index arrays in a cell array.
%       weight  - A double array of the weights of each of the variable
%                 above, indicating variables next to each other. Group 
%                 multiple SOS via cell arrays as above.
%
%   Quadratic Constraints (QC) [qrl <= x'Qx + l'x <= qru]:
%       Q       - A sparse double matrix of the quadratic terms for the
%                 constraint, group multiple quadratic constraints via a
%                 cell array of matrices. [NOT TRIL / TRIU]
%       l       - A column vector of the linear terms for the constraint,
%                 group multiple quadratic constraints in a matrix, each
%                 column representing each constraint.
%       qrl     - A scalar representing the quadratic constraint lower
%                 bound. Group multiple quadratic constraints in a
%                 column vector, each row representing each constraint.
%       qru     - A scalar representing the quadratic constraint upper 
%                 bound. Group multiple quadratic constraints in a
%                 column vector, each row representing each constraint.
%
%   Nonlinear Objective & Constraints (min f(x) and cl <= c(x) <= cu)
%       See opti_scipnl for usage. Contains multiple fields of instructions
%       lists which the MEX interface parses into SCIP expressions.
%
%   * From OPTI v1.78 scip can also solve AMPL .nl models. Simply pass the FULL path
%	  of the AMPL .nl file to scip as the first argument, and options (if requred) as
%	  the second argument. This interface will also write an AMPL solution file to the
%	  same directory as the input file. Remember SCIP does not solve models with 
%	  trigonometric functions.
%
%   Copyright (C) 2012-2013 Jonathan Currie (IPL)