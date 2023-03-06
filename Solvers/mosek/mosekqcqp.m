function [x,fval,exitflag,info] = mosekqcqp(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,x0,opts)
%MOSEKQCQP Solve a QCQP using MOSEK
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         qrl <= x'Qx + l'x <= qru
%                                              lb <= x <= ub
%                                              
%
%   x = mosekqcqp(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub) solves a QCQP where H and f 
%   are the objective matrix and vector respectively, A,rl,ru are the 
%   linear constraints, Q,l,qrl,qru are the quadratic constraints (see below) 
%   and lb,ub are the bounds.
%
%   x = mosekqcqp(H,...,ub,x0) uses x0 as the initial solution guess.
%
%   x = mosekqcqp(H,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver. 
%
%   [x,fval,exitflag,info] = mosekqcqp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   Quadratic Constraint Form:
%       - Single Quadratic Constraint
%           Q is matrix
%           l is a column vector
%           qrl is a scalar
%           qru is a scalar
%       - Multiple Quadratic Constraints
%           Q is a row cell array, each cell is a quadratic constraint Q
%           l is a matrix, each column for each quadratic constraint l
%           qrl is a column vector, each row for each quadratic constraint qrl
%           qru is a column vector, each row for each quadratic constraint qru

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

t = tic;

% Handle missing arguments
if nargin < 13, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 12, x0 = []; end
if nargin < 11, ub = []; end
if nargin < 10, lb = []; end
if nargin < 9, qru = []; end
if nargin < 8, qrl = []; end
if nargin < 7, l = []; end
if nargin < 6, Q = []; end
if nargin < 5, error('You must supply at least 5 arguments to mosekqcqp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
prob = mosekProb(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,[],[],x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);