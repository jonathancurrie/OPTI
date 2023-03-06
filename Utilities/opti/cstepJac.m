function jac = cstepJac(fun,x,nrow,h)
%CSTEPJAC  Complex Step Jacobian of a Scalar or Column Vector Function
%
%  cstepJac uses a complex step approach to finite difference to ensure
%  robust estimate of the Jacobian. Note your function MUST be
%  written in Matlab code only and must not contain any if/else statements.
%
%  If it does not work - ensure you are using .' (transpose) and not ' (complex conjugate).
%
%   jac = cstepJac(fun,x) calculates the Jacobian of fun at the point x.
%
%   jac = cstepJac(fun,x,nrow) specifies the number of rows in the output
%   of fun (i.e. fun is a vector function)
%
%   jac = cstepJac(fun,x,nrow,h) specifies the step-size. This defaults to
%   1e-8.
%
%
%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

if(nargin < 4), h = 1e-8; end
if(nargin < 3), nrow = []; end

if(numel(x) > 1)
    n = numel(x);
    if(isempty(nrow))
       [nrow,c] = size(fun(x));
       if(c > nrow), error('This function requires the output of fun to be a column vector'); end       
    end
    jac = zeros(nrow,n);
    for k=1:n
        xx = x; xx(k) = xx(k) + 1i*h;
        jac(:,k) = imag(fun(xx))./h;
    end
else
    jac = imag(fun(x + 1i*h))./h;
end
