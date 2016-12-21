function [x,fval,exitflag,info,Opt] = opti_lsqcurvefit(fun,x0,xdata,ydata,lb,ub,opts)
%OPTI_LSQCURVEFIT Solve a NLS using an OPTI NLS Solver (Matlab Overload)
%
%   [x,fval,exitflag,info] = opti_lsqcurvefit(fun,x0,xdata,ydata,lb,ub) 
%   solves the nonlinear least squares problem sum[(f(x,xdata)-ydata)^2] 
%   where fun is the nonlinear function to be fitted [fun(x,xdata)],
%   starting at x0, to the sample data ydata. Optional bounds lb and ub can
%   be placed on the decision variables x.
%
%   [x,fval,exitflag,info] = opti_lsqcurvefit(fun,...,ub,opts) allows the 
%   user to specify optiset options. This includes specifying a solver via 
%   the 'solver' field of optiset.
%
%   [x,...,info,Opt] = opti_lsqcurvefit(fun,...) returns the internally 
%   built OPTI object.

%   Copyright (C) 2011 Jonathan Currie (IPL)


% Handle missing arguments
if nargin < 7, opts = optiset; end 
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, error('You must supply at least 4 arguments to opti_lsqcurvefit'); end

%Build OPTI Object
Opt = opti('fun',fun,'data',xdata,ydata,'bounds',lb,ub,'x0',x0,'options',opts);

%Solve
[x,fval,exitflag,info] = solve(Opt);
