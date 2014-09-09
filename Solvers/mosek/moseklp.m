function [x,fval,exitflag,info] = moseklp(f,A,rl,ru,lb,ub,x0,opts)
%MOSEKLP Solve a LP using MOSEK
%
%   min f'*x      subject to:     rl <= A*x <= ru
%    x                            lb <= x <= ub                                
%
%   x = moseklp(f,A,b,Aeq,beq,lb,ub) solves a LP where f is the objective 
%   vector, A,rl,ru are the linear constraints and lb,ub are the bounds.
%
%   x = moseklp(f,...,ub,x0) uses x0 as the initial solution guess.
%
%   x = moseklp(f,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver. 
%
%   [x,fval,exitflag,info] = moseklp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 8, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 7, x0 = []; end
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, error('You must supply at least 4 arguments to moseklp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
prob = mosekProb([],f,A,rl,ru,[],[],[],[],lb,ub,[],[],x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);