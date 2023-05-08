% CBC  Solve a MILP/MIQP using CBC
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_cbc() INSTEAD!
%
% cbc uses the Coin-Or Branch and Cut library.
%
%   [x,fval,exitflag,nodes,xc] = cbc(H, f, A, rl, ru, lb, ub, xtype, sos, x0, opts)
%
%   Input arguments:
%       H - quadratic objective matrix (sparse, tril, optional) [NOT CURRENTLY SUPPORTED]
%       f - linear objective vector
%       A - linear constraint matrix (sparse)
%       rl - linear constraint lhs
%       ru - linear constraint rhs
%       lb - decision variable lower bounds (optional)
%       ub - decision variable upper bounds (optional)
%       xtype - decision variable integrality ('C', 'I' or 'B')
%       sos - SOS structure with fields type, index and weight
%       x0 - initial solution guess
%       opts - solver options (see below)
%
%   Return arguments:
%       x - integer solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       nodes - number of nodes searched by the solver
%       xc - continuous solution vector
%
%   Option Fields (all optional - see cbcset):
%       intTol - absolute integer tolerance
%       maxnodes - maximum nodes to explore
%       maxtime - maximum execution time [s]
%       display - solver display level [0,1]
%
%   Return Status:
%       0 - looks optimal
%       1 - lp relaxation infeasible
%       2 - gap reached
%       3 - maximum nodes reached
%       4 - maximum time reached
%       5 - user exited
%       6 - #solutions reached
%       7 - lp relaxation unbounded
%		8 - proven infeasible
%
%
%   Copyright (C) 2013 Jonathan Currie (Control Engineering)
%
%   See Also opti_cbc