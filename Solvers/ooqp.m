% OOQP  Solve a LP or QP using OOQP
%
% THIS IS A LOW LEVEL FUNCTION - USE opti_ooqp() INSTEAD!
%
% ooqp uses the Object Orientated Quadratic Programming library.
%
%   [x,fval,stat,iter,lambda] = ooqp(H, f, A, rl, ru, Aeq, beq, lb, ub, opts)
%
%   Input arguments:
%       H - quadratic objective matrix (sparse, tril, optional)
%       f - linear objective vector
%       A - linear inequality matrix 
%       rl - linear inequality lower bounds
%       ru - linear inequality upper bounds
%       Aeq - linear equality matrix
%       beq - linear equality rhs
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       opts - solver options (see below)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       stat - exit status (see below)
%       iter - number of iterations taken by the solver
%       lambda - structure of dual information
%
%   Option Fields (all optional - see ooqpset):
%       maxiter - maximum solver iterations
%       maxtime - maximum solver execution time
%       display - solver display level [0,1,2]
%       objbias - Objective bias term (defaults to 0.0)
%       algorithm - Solver algorithm (Mehrotra 0, Gondzio {1})
%       linear_solver - Linear solver algorithm (PARDISO 0, MA57 {1}, MA27 2)
%
%   Return Status:
%       0 - optimal
%       1 - not finished
%       2 - maximum iterations exceeded
%       3 - infeasible
%       4 - ooqp error
%
%
%   Code is based in parts on original MEX interface by E. Michael Gertz, 
%   Stephen J. Wright
%
%   Copyright (C) 2013 Jonathan Currie (Control Engineering)
%
%   See Also opti_ooqp ooqpset