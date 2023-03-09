function [x,fval,exitflag,info] = opti_levmar(fun,grad,x0,ydata,lb,ub,A,b,Aeq,beq,opts)
%OPTI_LEVMAR Solve a NLS using LEVMAR (Levenberg-Marquardt by Manolis Lourakis)
%
%   min sum[ (F(x) - ydata)^2 ]       subject to:   A*x <= b
%    x                                              Aeq*x = beq
%                                                   lb <= x <= ub
%
%   x = opti_levmar(fun,grad,x0,ydata) solves a NLS where fun is the 
%   fitting function. grad is an optional gradient of the fitting function 
%   and x0 is a starting guess. ydata is the data to fit the function to. 
%
%   x = opti_levmar(fun,grad,x0,ydata,lb,ub) solves subject to decision
%   variables bounds lb <= x <= ub. Avoid Infinite bounds.
%
%   x = opti_levmar(fun,...,ub,A,b) solves subject to the linear
%   inequalities Ax <= b.
%
%   x = opti_levmar(fun,...,b,Aeq,beq) solves subject to the linear
%   equalities Aeqx = beq.
%
%   x = opti_levmar(fun,...,beq,opts) uses opts to pass optiset options to 
%   the solver. 
%
%   [x,fval,exitflag,info] = opti_levmar(...) returns the objective value 
%   at the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR LEVMAR
%   See referenced GNU Public License

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

if(nargin < 11), opts = optiset; end
if(nargin < 10), beq = []; end
if(nargin < 9), Aeq = []; end
if(nargin < 8), b = []; end
if(nargin < 7), A = []; end
if(nargin < 6), ub = []; end
if(nargin < 5), lb = []; end
if(nargin < 4), error('LEVMAR requires at least 4 arguments'); end

%Setup display level
opts.display = dispLevel(opts.display);
opts.optiver = optiver;

%Check we have a valid x0
if(isempty(x0) || any(isnan(x0)))
    error('LEVMAR requires an initial guess, x0!');
end

t = tic;
% Run LEVMAR
[x, fval, solverflag, iter, feval] = levmar(fun,grad,x0,ydata,lb,ub,A',b,Aeq',beq,opts);

%Collect Results
info.Iterations = iter;
info.FuncEvals = feval;
info.Time = toc(t);
info.Algorithm = 'LEVMAR: Levenberg-Marquardt in C/C++';

switch(solverflag)
    case {1,2}
        info.Status = 'Optimal'; exitflag = 1;
    case 3
        info.Status = 'Exceeded Iterations'; exitflag = 0;
    case 4
        info.Status = 'Error - Singular Matrix Detected'; exitflag = -2;
    case 5
        info.Status = 'Error - No Further Reduction Possible'; exitflag = -2;
    case 6
        info.Status = 'Error - Stopped By Small ||e||_2'; exitflag = -2;
    case 7
        info.Status = 'Error - NaN or Inf Detected'; exitflag = -3;
    otherwise        
        info.Status = 'LEVMAR Error'; exitflag = -3;
end
