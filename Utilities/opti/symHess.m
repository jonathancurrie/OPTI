function [hess,pattern,symhess] = symHess(fun,nvar,isTril,var)
%SYMHESS  Returns a Symbolically Differentiated Hessian
%
%   hess = symHess(fun) uses the symbolic toolbox to automatically generate
%   the Hessian of the function handle fun. You should supply the function
%   handle in the form @(x) x(1)^2 + x(2), noting x must be the only
%   variable used, and is indexed in the equation (no vector operations).
%
%   hess = symHess(fun,nvar) specifies the number of variables in the
%   equation, assuming consecutive ordering. Useful if a variable is not
%   specified in the original equation to pad the Hessian with zeros.
%
%   hess = symHess(fun,nvar,isTril) specifies if the returned Hessian
%   should be Symmetric Lower Triangular.
%
%   [hess,pattern] = symHess(...) also returns the sparsity pattern as a
%   function handle.

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 4), var = 'x'; end
if(nargin < 3), isTril = false; end
if(nargin < 2 || isempty(nvar)), nvar = 0; end

if(~optiCheckSymTBX())
    hess = []; symhess = [];
    return
end

if(~isa(fun,'function_handle') && ~isa(fun,'barvec'))
    error('Fun should be a function handle!');
end
if(isa(fun,'function_handle') && nargin(fun) ~= 1)
    error('fun should only have one input argument (x)');
end

%Convert function to a symbolic expression
if(isa(fun,'function_handle'))
    [symfun,ind] = func2sym(fun,var);
    ind = ind{1};
else
    symfun = sym(getEq(fun)); %convert from barvec to symbolic expression
    v = char(symvar(symfun));
    ind = 1:length(strfind(v,var));
end

%Check Symbolic Expression is Scalar
if(~isscalar(symfun))
    error('symHess only works on scalar functions');
end

%Create each symbolic partial derivative
symhess = symPartialDer(symfun,var,nvar,ind,true);

%Convert if required to lower triangular
if(isTril), symhess = tril(symhess); end

%Extract Sparsity Pattern
sp = zeros(size(symhess)); sp(symhess ~= 0) = 1; sp = sparse(sp);
pattern = @() sp;

%Return to function handle
hess = sym2func(symhess,var);
