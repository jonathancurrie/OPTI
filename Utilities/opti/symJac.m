function [jac,pattern,symjac] = symJac(fun,nvar,tran,var,file)
%SYMJAC  Returns a Symbolically Differentiated Jacobian
%
%   jac = symJac(fun) uses the symbolic toolbox to automatically generate
%   the Jacobian of the function handle fun. You should supply the function
%   handle in the form @(x) x(1)^2 + x(2), noting x must be the only
%   variable used, and is indexed in the equation (no vector operations).
%
%   jac = symJac(fun,nvar) specifies the number of variables in the
%   equation, assuming consecutive ordering. Useful if a variable is not
%   specified in the original equation to pad the Jacobian with zeros.
%
%   jac = symJac(fun,nvar,tran) specifies to transpose (.') the generated
%   Jacobian prior to returning the function handle.
%
%   jac = symJac(fun,nvar,tran,var) specifies the variable to look for as a
%   string. i.e. instead of x(1)^2, you could use z(1)^2 and supply 'z'.
%
%   jac = symJac(fun,nvar,tran,var,file) specifies a file to write the
%   jacobian to. A .m extension specifies a MATLAB file, a .c extension
%   specifies a c file (TBC).
%
%   [jac,pattern] = symJac(...) also returns the sparsity pattern as a
%   function handle.

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 5), file = []; end
if(nargin < 4 || isempty(var)), var = 'x'; end    
if(nargin < 3 || isempty(tran)), tran = 0; end
if(nargin < 2 || isempty(nvar)), nvar = 0; end

if(~optiCheckSymTBX())
    jac = []; pattern = []; symjac = [];
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

%Create each symbolic partial derivative
symjac = symPartialDer(symfun,var,nvar,ind);

%Transpose if required
if(tran), symjac = symjac.'; end

%If we have a file, write to it and return a handle to the file
if(~isempty(file))
    %Decide if .m or .c
    if(~isempty(strfind(file,'.m')))
        jac = matlabFunction(symjac,'file',file,'vars',{symvar(symjac).'});
    elseif(~isempty(strfind(file,'.c')))
        error('not yet implemented');
        jac = symMakeCMex(symjac,file);        
    end
else
    %Return to function handle
    jac = sym2func(symjac,var);
end

%Extract Sparsity Pattern
sp = zeros(size(symjac)); sp(symjac ~= 0) = 1; sp = sparse(sp);
pattern = @() sp;





