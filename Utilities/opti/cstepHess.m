function hess = cstepHess(fun,x,isTril,h)
%CSTEPHESS  Complex Step Hessian of a Scalar Function
%
%  cstepHess uses a complex step approach to finite difference to ensure
%  robust estimate of the Hessian. Note your function MUST be
%  written in Matlab code only and must not contain any if/else statements.
%
%  If it does not work - ensure you are using .' (transpose) and not ' (complex conjugate).
%
%   hess = cstepHess(fun,x) calculates the Hessian of fun at the point x.
%
%   hess = cstepHess(fun,x,isTril) specifies if the returned Hessian
%   should be Symmetric Lower Triangular.
%
%   hess = cstepHess(fun,x,isTril,h) specifies the step-size. This defaults to
%   1e-3.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)
%
%   This code follows ideas from "New Complex-Step Derivative Approximations
%   with Application to Second-Order Kalman Filtering", by Kok-Lam Lai, John
%   L. Crassidis and Yang Cheng.

if(nargin < 4), h = 1e-3; end
if(nargin < 3 || isempty(isTril)), isTril = false; end

%Constants from paper
I = sqrt(2)/2*(1i + 1);
n = length(x); hess = zeros(n,n);
%Diagonal Elements
for k = 1:n
    xu = x; xu(k) = xu(k) + I*h;
    xl = x; xl(k) = xl(k) - I*h;
    hess(k,k) = imag((fun(xu) + fun(xl))/h^2); 
end
%Lower Triangular Elements
lam = 1; k = n-1;
while(k > 0)
    for phi = 1:k
        xu = x; xu(phi:phi+lam) = xu(phi:phi+lam) + I*h;
        xl = x; xl(phi:phi+lam) = xl(phi:phi+lam) - I*h;
        Fsum = sum(sum(hess(phi:phi+lam,phi:phi+lam)));
        hess(phi+lam,phi) = (imag((fun(xu) + fun(xl))/h^2) - Fsum)/2;
        hess(phi,phi+lam) = hess(phi+lam,phi); %copy to upper tri        
    end
    k = k - 1;
    lam = lam + 1;
end
%Make tril if requested
if(isTril), hess = tril(hess); end