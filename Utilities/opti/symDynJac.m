function [dfdz,dfdp,sym_dfdz,sym_dfdp] = symDynJac(ode,nstates,nparam)
%SYMJAC  Returns Symbolically Differentiated Partial Derivatives of a Dynamic System
%
%   [dfdz,dfdp] = symDynJac(ode) uses the symbolic toolbox to automatically 
%   generate the sensitivity partial derivatives of the function handle ode. 
%   You should supply the function handle in the form @(t,z,p) z(1)^2 + p(1)
%   noting t, z and p must be the only variables used, and are indexed in the 
%   equation (no vector operations).
%
%   [dfdz,dfdp] = symDynJac(ode,nstates) specifies the number of states in 
%   the equation, assuming consecutive ordering. Useful if a state is not
%   specified in the original equation to pad DFDZ with zeros.
%
%   [dfdz,dfdp] = symDynJac(ode,nstates,nparam) specifies the number of 
%   parameters in the equation, assuming consecutive ordering. Useful if a 
%   parameter is not specified in the original equation to pad DFDP with 
%   zeros.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

if(nargin < 3), nparam = 0; end
if(nargin < 2), nstates = 0; end

if(~optiCheckSymTBX())
    dfdz = []; sym_dfdz = [];
    dfdp = []; sym_dfdp = [];
    return
end

if(~isa(ode,'function_handle') && ~isa(ode,'barvec'))
    error('Fun should be a function handle!');
end
if(isa(ode,'function_handle') && nargin(ode) ~= 3)
    error('ODE should only have three input arguments (t,z,p)');
end

%Convert ODE to a symbolic expression
if(isa(ode,'function_handle'))
    [symode,ind] = func2sym(ode,{'z','p','t'});
    indz = ind{1}; indp = ind{2};
else
    symode = sym(getEq(ode)); %convert from barvec to symbolic expression
    v = char(symvar(symode));
    indz = 1:length(strfind(v,'z'));
    indp = 1:length(strfind(v,'p'));
end

%Create each symbolic partial derivative
sym_dfdz = symPartialDer(symode,'z',nstates,indz);
sym_dfdp = symPartialDer(symode,'p',nparam,indp);

%Return to function handles (note order of args important and changed)
dfdz = sym2func(sym_dfdz,{'t','z','p'});
dfdp = sym2func(sym_dfdp,{'t','z','p'});


