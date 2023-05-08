function [x,fval,exitflag,info] = moseksdp(f,A,rl,ru,lb,ub,sdcone,x0,opts)
%MOSEKSDP Solve a SDP using MOSEK
%
%   min f'*x      subject to:     rl <= A*x <= ru
%    x                            lb <= x <= ub
%                                 X = F1*x1 + F2*x2 + ... + Fn*xn - F0
%                                 X >= 0 [positive semidefinite]
%                                              
%   x = moseksdp(f,A,rl,ru,lb,ub,sdcone) solves a SDP where 
%   f is the objective vector, A,rl,ru are the linear constraints, lb,ub 
%   are the bounds and sdcone contains the semidefinite constraints. Each 
%   cell within sdcone is a sparse matrix of the form [F0(:) F1(:) F2(:)...] 
%   (or in alternative notation, [C A0 A1 A2..]) where each F matrix has 
%   been converted to a column vector and concatenated.
%
%   x = moseksdp(f,...,sdcone,x0) uses x0 as the initial solution guess.
%
%   x = moseksdp(H,...,x0,opts) uses opts to pass mosekset options to the
%   MOSEK solver.
%
%   [x,fval,exitflag,info] = moseksdp(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

t = tic;

% Handle missing arguments
if nargin < 9, opts = mosekset('warnings','off'); else opts = mosekset(opts); end 
if nargin < 8, x0 = []; end
if nargin < 7, sdcone = []; end
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, error('You must supply at least 4 arguments to moseksdp'); end

%Build MOSEK Command
[cmd,param] = mosekBuild(opts);

%Create Problem Structure
prob = mosekProb([],f,A,rl,ru,[],[],[],[],lb,ub,[],sdcone,x0,opts);
    
%Call MOSEK
[rcode,res] = mosekopt(cmd,prob,param);

%Extract Results
[x,fval,exitflag,info] = mosekRes(prob,rcode,res,t);