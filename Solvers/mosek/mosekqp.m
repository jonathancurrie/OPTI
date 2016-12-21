function [x,fval,exitflag,info] = mosekqp(H,f,A,rl,ru,lb,ub,x0,opts)
%MOSEKQP Solve a QP using MOSEK
%
%   min 0.5*x'*H*x + f'*x      subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                              
%   x = mosekqp(H,f,A,b,Aeq,beq,lb,ub) solves a QP where H and f are the 
%   objective matrix and vector respectively, A,b are the linear constraints
%   and lb,ub are the bounds.
%
%   x = mosekqp(H,...,ub,x0) uses x0 as the initial solution guess.
%
%   x = mosekqp(H,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver.
%
%   [x,fval,exitflag,info] = mosekqp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 9, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 8, x0 = []; end
if nargin < 7, ub = []; end
if nargin < 6, lb = []; end
if nargin < 5, error('You must supply at least 5 arguments to mosekqp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
prob = mosekProb(H,f,A,rl,ru,[],[],[],[],lb,ub,[],[],x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);