% ASL  Read an AMPL .NL file and convert to Matlab Values
%
%   NOTE - This API has substantially changed from OPTI v1.79
%
% THIS IS A LOW LEVEL FUNCTION - USE amplRead() INSTEAD!
%
% asl uses the Netlib AMPL Solver Library to open and parse the file. The
% MEX file is a little different in that it stays open after a call to
% asl(), in order to evaluate nonlinear functions later. You must use the
% close command in order to close the mex file.
%
%   [varargout] = asl(command,varargin)
%
%   Input arguments:
%       command - a string containing the command to execute (see below)
%       varargin - command specific arguments
%
%   Return arguments:
%       varagout - command specific return variables
%
%   Command: 'open' [Open an AMPL file and read basic information + check if LP/QP/QCQP]
%       [aslprob,sizes] = asl('open',path,isNLP)
%
%       path - full path to the file on your PC
%       isNLP - 1 to skip identification of LP/QP/QCQP (assume (MI)NLP)
%
%       aslprob structure fields:
%           H  - Quadratic Objetive H matrix (sparse, if QP)
%           f  - Linear Objective f vector (if LP/QP)
%           lb - Decision variable lower bounds
%           ub - Decision variable upper bounds    
%           A  - Linear Constraints A matrix (sparse, if LP/QP)
%           cl - Constraints lower bounds
%           cu - Constraints upper bounds
%           Q  - Cell array of quadratic constraint Q matrices (sparse, if QCQP)
%           l  - Matrix of quadratic constraint linear vectors (if QCQP)
%           qcind - Vector of indices of Quadratic Constraints (to remove from A and cl,cu, if QCQP)
%           x0 - initial decision variable guess
%           v0 - initial dual variable guess
%           sense - 1: minimization, -1: maximization
%           conlin - vector of constraint linearity (0 linear,<0  nonlinear, >0 quadratic)           
%
%       sizes - vector of problem sizes (see amplRead.m)
%
%   Command: 'isopen' [Check if AMPL ASL interface has been opened]
%       asl('isopen')
%
%   Command: 'close' [Close an opened AMPL file and destory MEX file memory]
%       asl('close')
%
%   Command: 'fun' [Evaluate objective function at x]
%       f = asl('fun',x)
%
%   Command: 'grad' [Evaluate objective gradient at x]
%       g = asl('grad',x)
%
%   Command: 'con' [Evaluate constraint function at x]
%       c = asl('con',x)
%
%   Command: 'jac' [Evaluate constraint Jacobian at x]
%       j = asl('jac',x)
%       js = asl('jac',x,1) %return sparse Jacobian
%
%   Command: 'jacstr' [Evaluate constraint Jacobian Sparsity Structure]
%       s = asl('jacstr')
%
%   Command: 'hess' [Evaluate Hessian of the Lagrangian (L = sigma*Hobj + sum(v*Hcon))]
%       h = asl('hess',x,sigma,lambda)
%       hs = asl('hess',x,sigma,lambda,1) %return sparse Hessian
%
%   Command: 'hessstr' [Evaluate Hessian Sparsity Structure]
%       s = asl('hessstr')
%
%   Copyright (C) 2011-2013 Jonathan Currie (I2C2)