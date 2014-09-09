% NOMAD  Solve a Global MINLP/NLP using NOMAD, a Blackbox Optimization Library
%
%   min f(x)      subject to:     nlcon(x) <= nlrhs
%    x                            lb <= x <= ub
%                                 xi in Z
%                                 xj in {0,1}
%            
%
%   [x,fval,exitflag,iter,nfval] = nomad(fun,x0,lb,ub,nlcon,nlrhs,xtype,opts)
%
%   Input arguments:
%       fun - nonlinear function handle 
%       x0 - initial solution guess
%       lb - decision variable lower bounds
%       ub - decision variable upper bounds
%       nlcon - nonlinear inequalities function handle
%       nlrhs - nonlinear inequalities rhs 
%       xtype - integer variable string ('C' - continuous, 'I' - integer, 'B', binary)
%       opts - solver options (see nomadset)
%
%   Return arguments:
%       x - solution vector
%       fval - objective value at the solution
%       exitflag - exit status (see below)
%       iter - number of iterations taken by the solver
%       nfval - number of function evaluations taken by the solver
%
%   Return Status:
%       1 - converged / target reached
%       0 - maximum iterations / function evaluations exceeded
%      -1 - infeasible / mesh limit reached
%      -2 - initialization error
%      -3 - nomad error
%      -5 - user exit
%
%
%   NOMAD is released under the Lesser GNU Public License (LGPL).
%   Type nomad('-info') to see license and author details.
%
%   MEX Interface Copyright (C) 2012 Jonathan Currie (I2C2)