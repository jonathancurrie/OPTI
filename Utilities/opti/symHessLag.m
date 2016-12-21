function [hess,pattern,symhess] = symHessLag(obj,nlcon,nvar,isTril,var)
%SYMHESSLAG  Returns a Symbolically Differentiated Hessian of the Lagrangian
%
%   hess = symHessLag(obj,nlcon) uses the symbolic toolbox to automatically 
%   generate the Hessian of the Lagrangian of the objective function obj and
%   nonlinear constraints nlcon. You should supply the function
%   handles in the form @(x) x(1)^2 + x(2), noting x must be the only
%   variable used, and is indexed in the equation (no vector operations).
%
%   hess = symHessLag(obj,nlcon,nvar) specifies the number of variables in 
%   the equation, assuming consecutive ordering. Useful if a variable is not
%   specified in the original equation to pad the Hessian with zeros.
%
%   hess = symHessLag(obj,nlcon,nvar,isTril) specifies if the returned 
%   Hessian should be Symmetric Lower Triangular.
%
%   [hess,pattern] = symHessLag(...) also returns the sparsity pattern as a
%   function handle.

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 5), var = 'x'; end
if(nargin < 4), isTril = false; end
if(nargin < 3 || isempty(nvar)), nvar = 0; end
if(nargin < 2), nlcon = []; end

if(~optiCheckSymTBX())
    hess = []; symhess = [];
    return
end

if(~isa(obj,'function_handle') && ~isa(obj,'barvec'))
    error('Objective should be a function handle!');
end
if(isa(obj,'function_handle') && nargin(obj) ~= 1)
    error('Objective should only have one input argument (x)');
end
%If no constraints, then easy one - just use symHess
if(isempty(nlcon))
    [~,symhess] = symHess(obj,nvar,isTril,var);
    %Multiply by sigma
    symhess = sym('sigma')*symhess;
    %Return to function handle with 3 args
    hess = sym2func(symhess,sprintf('%s,sigma,lambda',var));
    return;
end
%Otherwise, standard problem
if(~isa(nlcon,'function_handle') && ~isa(nlcon,'barvec'))
    error('Nonlinear Constraints should be a function handle!');
end
if(isa(nlcon,'function_handle') && nargin(nlcon) ~= 1)
    error('Nonlinear Constraints should only have one input argument (x)');
end

%Convert objective to a symbolic expression
if(isa(obj,'function_handle'))
    [symobj,ind] = func2sym(obj,var);
    indo = ind{1};
else
    symobj = sym(getEq(obj)); %convert from barvec to symbolic expression
    v = char(symvar(symobj));
    indo = 1:length(strfind(v,var));
end

%Check Symbolic Objective is Scalar
if(~isscalar(symobj))
    error('The Objective Function must return a scalar');
end

%Convert constraints to a symbolic expression
if(isa(nlcon,'function_handle'))
    symcon = func2sym(nlcon,var);
else
    symcon = sym(getEq(nlcon)); %convert from barvec to symbolic expression
end

%If we don't know number of vars, then use max of what is found so far
if(nvar == 0)
    nvar = max(numel(symvar(symobj)),numel(symvar(symcon)));
end

%Begin creating hessLag
symhess = sym('sigma')*symPartialDer(symobj,var,nvar,indo,true);

%For each constraint equation, multiply by lambda(i) then add to our hessian
for i = 1:length(symcon)
    symhess = symhess + sym(sprintf('lambda(%d)',i))*symPartialDer(symcon(i),var,nvar,[],true);
end

%Convert if required to lower triangular
if(isTril), symhess = tril(symhess); end

%Extract Sparsity Pattern
sp = zeros(size(symhess)); sp(symhess ~= 0) = 1; sp = sparse(sp);
pattern = @() sp;

%Return to function handle
hess = sym2func(symhess,sprintf('%s,sigma,lambda',var));
