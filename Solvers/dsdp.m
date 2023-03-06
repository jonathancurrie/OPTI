% DSDP  Solve a SDP using DSDP [Alternative Dual Form]
%
% SDP Form:
%  (P)  MINIMIZE trace(C*X) SUCH THAT  trace(A_i*X) = f_i, i=1,...,m and X >= 0
%  (D)  MAXIMIZE   dot(f,y) SUCH THAT  C - SUM A_i*y_i >= 0.
%   
%   To supply in standard primal form simply negate f and each matrix in 
%   each semidefinite cone.
%
%
% dsdp uses the DSDP library, (C) University of Chicago.
%
%   [y,fvals,exitflag,stats,X] = dsdp(f,A,b,lb,ub,sdcone,y0,opts)
%
%   Input arguments:
%       f - linear objective vector 
%       A - linear constraint matrix (sparse)
%       b - linear constraint rhs (A*y <= b)
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds (lb <= y <= ub)
%       sdcone - Semidefinite Cone constraints (see below)
%       y0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       y - solution vector
%       fval - structure with pval (primal) and dval (dual) objective values at the solution
%       exitflag - exit status (see below)
%       stats - solution status structure
%       X - cell array of solution vectors of the form [LP;SDP...;FIX] (use opti_dmat to reform SDP X to matrices)
%
%   Semidefinte Cones:
%       sdcone is a cell array where each cell is a semidefinite cone. Each 
%       cone is a sparse matrix of the form:
%           sdcone{1} = [C A0 A1 A2]
%       where each matrix is assumed to be square and symmetric, thus only
%       the upper triangular elements are supplied as a column each. For 
%       example, if C and A0 are 2x2 symmetric matrices in a problem with 
%       one variable, then supplied to dsdp would be:
%           ind = triu(ones(2))==1;
%           sdcone{1} = sparse([C(ind) A0(ind)]);
%
%   options Fields:
%       maxiter - maximum solver iterations
%       maxtime - maximum solver execution time
%       display - 0 to 2 display level
%       For the remainder see dsdpset.
%
%   Return Status:
%       1 - dsdp converged
%       5 - dual objective big enough to stop
%       7 - user termination
%      -2 - numerical difficulties with short step lengths
%      -3 - maximum number of iterations
%      -6 - initial points imply S is not positive
%      -8 - indefinite schur matrix
%      -9 - numerical error
%     -27 - maximum time reached
%
%   Also Check stats.pdflag:
%       0 - not sure if primal or dual feasible, check bounds
%       1 - both primal and dual feasible and bounded
%       3 - dual unbounded, primal infeasible
%       4 - dual infeasible, primal unbounded
 
%   Based in parts on the original MEX interface by DSDP authors.
%
%   Copyright (C) 2013 Jonathan Currie (Control Engineering)