% CLP  Solve a LP or QP using CLP
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_clp() INSTEAD!
%
% clp uses the Coin-Or Linear Programming library.
%
%   [x,fval,exitflag,iter,lambda] = clp(H, f, A, rl, ru, lb, ub, opts)
%
%   Input arguments:
%       H - quadratic objective matrix (sparse, tril, optional)
%       f - linear objective vector
%       A - linear constraint matrix (sparse)
%       rl - linear constraint lhs
%       ru - linear constraint rhs
%       lb - decision variable lower bounds (optional)
%       ub - decision variable upper bounds (optional)
%       opts - solver options (see below)      
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       lambda - structure of dual information
%
%   Option Fields (all optional - see clpset):
%       algorithm       - Solver algorithm (DualSimplex 0, PrimalSimplex 1, PrimalSimplexOrSprint 2, Barrier 3, BarrierNoCross 4, Automatic {5})
%       primalTol       - Primal Tolerance
%       dualTol         - Dual Tolerance
%       maxiter         - maximum solver iterations
%       maxtime         - maximum execution time [s]
%       display         - solver display level [0,1->100 increasing]
%       objbias         - Objective bias term (defaults to 0.0)
%       numPresolvePasses - Number of presolver passes {5}
%       factorFreq      - Simplex factorization frequency (empty uses internal heuristic)
%       numberRefinements - Number of iterative Simplex refinements {0}
%       primalObjLim    - Primal objective limit {Inf}
%       dualObjLim      - Dual objective limit {Inf}
%       numThreads      - Number of Cilk Worker Threads (LP Only and > 1 only with Aboca CLP Build)
%       abcState        - Aboca Partition Size (empty uses internal heuristic)
%
%   Return Status:
%       0 - Proven Optimal
%       1 - Proven Primal Infeasible
%       2 - Proven Dual Infeasible
%       3 - Maximum Iterations or Time Reached
%       4 - Errors Present
%       5 - User Exited
%
%
%   Copyright (C) 2013 Jonathan Currie (I2C2)
%
%   See Also opti_clp