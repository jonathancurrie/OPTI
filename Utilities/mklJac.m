% MKLJAC  Estimate the Jacobian of a function via Numerical Differences
%
% mklJac uses the Intel Math Kernel Library (MKL) djacobi function to
% estimate the Jacobian of a function using central differences.
%
%   jac = mklJac(fun,x) uses the supplied function handle fun and the
%   current state vector x to estimate the gradient of the function.
%
%   jac = mklJac(fun,x,nrow) specifies the number of rows in the vector
%   returned from fun, used to determine the return size of the Jacobian.
%   If the number of rows is not specified, then a dummy function call is
%   used in order to determine the number of rows.
%
%   jac = mklJac(fun,x,nrow,tol) specifies the tolerance of the numerical
%   difference algorithm of all variables. By default tol is 1e-6.
%
%   [jac,status] = mklJac(fun,x,nrow,tol) also returns 1 if the algorithm
%   was successful, or 0 if it failed.
%
%
%   Copyright (C) 2011 Jonathan Currie (I2C2)