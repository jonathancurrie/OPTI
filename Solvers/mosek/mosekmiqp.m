function [x,fval,exitflag,info] = mosekmiqp(H,f,A,rl,ru,lb,ub,xint,x0,opts)
%MOSEKMIQP Solve a MIQP using MOSEK
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1}                                              
%
%   x = mosekmiqp(H,f,A,rl,ru,lb,ub,xint) solves a QP where H and f are 
%   the objective matrix and vector respectively, A,rl,ru are the linear 
%   constraints, lb,ub are the bounds and xint is a string of integer 
%   variables ('C', 'I', 'B')
%
%   x = mosekmiqp(H,...,xint,x0) uses x0 as the initial solution guess.
%
%   x = mosekmiqp(H,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver.
%
%   [x,fval,exitflag,info] = mosekmiqp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 10, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 9, x0 = []; end
if nargin < 8, xint = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to mosekmiqp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
prob = mosekProb(H,f,A,rl,ru,[],[],[],[],lb,ub,xint,[],x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);