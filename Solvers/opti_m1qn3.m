function [x,fval,exitflag,info] = opti_m1qn3(fun,grad,x0,opts)
%OPTI_M1QN3 Solve a UNO using M1QN3 (Large-Scale L-BFGS Minimization)
%
%   min f(x)
%    x
%
%   x = opti_m1qn3(fun,grad,x0) solves a UNO where fun is the objective
%   function and grad is the gradient of the objective function. x0 
%   is a starting guess.
%
%   x = opti_m1qn3(fun,grad,x0,opts) uses opts to pass optiset 
%   options to the solver. 
%
%   [x,fval,exitflag,info] = opti_m1qn3(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR M1QN3
%   See supplied GPL v3 License

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(nargin < 4), opts = optiset; end
if(nargin < 3), error('M1QN3 requires at least 3 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('M1QN3 requires an initial guess, x0!');
end

t = tic;
% Run M1QN3
[x, fval, exitflag, iter, feval] = m1qn3(fun,grad,x0,opts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'M1QN3: Large-Scale L-BFGS Minimization';

switch(exitflag)
    case 1
        info.Status = 'Optimal';
    case {4,5,8}
        info.Status = 'Exceeded Iterations / Function Evaluations / Time';
        exitflag = 0;
    case {2,3,6,7}
        info.Status = 'M1QN3 Error';
        exitflag = -2;
    case 0
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise        
        info.Status = 'M1QN3 Error';
        exitflag = -3;
end
