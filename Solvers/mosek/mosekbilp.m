function [x,fval,exitflag,info] = mosekbilp(f,A,rl,ru,x0,opts)
%MOSEKBILP Solve a BILP using MOSEK
%
%   min f'*x      subject to:     rl <= A*x <= ru
%    x                            x in {0,1}
%                                 
%   x = mosekbilp(f,A,rl,ru) solves a BILP where f is the objective 
%   vector and A,rl,ru are the linear constraints
%
%   x = mosekbilp(f,...,ru,x0) uses x0 as the initial solution guess.
%
%   x = mosekbilp(f,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver.
%
%   [x,fval,exitflag,info] = mosekbilp(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (IPL)

t = tic;

% Handle missing arguments
if nargin < 6, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 5, x0 = []; end
if nargin < 4, error('You must supply at least 4 arguments to mosekbilp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
xint = repmat('B',size(f)); %all binary variables
prob = mosekProb([],f,A,rl,ru,[],[],[],[],[],[],xint,[],x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);