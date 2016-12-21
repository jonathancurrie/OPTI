function hess = cstepHessLag(obj,nlcon,x,sigma,lambda,isTril,h)
%CSTEPHESSLAG  Complex Step Hessian of the Lagrangian
%
%  cstepHessLag uses a complex step approach to finite difference to ensure
%  robust estimate of the Hessian of the Lagrangian. Note your function MUST be
%  written in Matlab code only and must not contain any if/else statements.
%
%  If it does not work - ensure you are using .' (transpose) and not ' (complex conjugate).
%
%   hess = cstepHessLag(obj,nlcon,x,sigma,lambda) uses complex step 
%   differentiation to automatically generate the Hessian of the Lagrangian 
%   of the objective function obj and nonlinear constraints nlcon.
%
%   hess = cstepHessLag(obj,nlcon,x,sigma,lambda,isTril) specifies if the 
%   returned Hessian should be Symmetric Lower Triangular.
%
%   hess = cstepHessLah(obj,...,isTril,h) specifies the step-size. This 
%   defaults to 1e-3.

%   Copyright (C) 2013 Jonathan Currie (IPL)
%
%   This code follows ideas from "New Complex-Step Derivative Approximations
%   with Application to Second-Order Kalman Filtering", by Kok-Lam Lai, John
%   L. Crassidis and Yang Cheng.

if(nargin < 7), h = 1e-3; end
if(nargin < 6 || isempty(isTril)), isTril = false; end

if(~isa(obj,'function_handle'))
    error('Objective should be a function handle!');
end
if(isa(obj,'function_handle') && nargin(obj) ~= 1)
    error('Objective should only have one input argument (x)');
end
%If no constraints, then easy one - just use cstepHess
if(isempty(nlcon))
    hess = sigma.*cstepHess(obj,x,isTril,h);
    return;
end
%Otherwise, standard problem
if(~isa(nlcon,'function_handle'))
    error('Nonlinear Constraints should be a function handle!');
end
if(isa(nlcon,'function_handle') && nargin(nlcon) ~= 1)
    error('Nonlinear Constraints should only have one input argument (x)');
end

%Begin creating hessLag
hess = sigma.*cstepHess(obj,x,isTril,h);

%For each constraint equation, multiply by lambda(i) then add to our hessian
parfor i = 1:length(lambda)
    hess = hess + lambda(i).*cstepHess(@(x) idxfnc(nlcon,x,i),x,isTril,h);
end

function g = idxfnc(fun,x,i)
g = fun(x); g = g(i);
