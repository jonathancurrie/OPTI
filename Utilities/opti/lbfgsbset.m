function options = lbfgsb(varargin)
%LBFGSBSET  Create or alter the options for Optimization with L-BFGS-B
%
% options = lbfgsbset('param1',value1,'param2',value2,...) creates an L-BFGS-B
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the L-BFGS-B
% default.
%
% options = lbfgsbset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = lbfgsbset() creates an options structure with all fields set to
% lbfgsbset defaults.
%
% lbfgsbset() prints a list of all possible fields and their function.
%
% See supplied L-BFGS-B Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'nupdate','pgtol'};
Defaults = {5,1e-5};        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)   
    %Scalar non negative double
    case {'pgtol'}
        err = opticheckval.checkScalarNonNeg(value,field);    
    %Scalar non negative integer
    case {'nupdate'}
        err = opticheckval.checkScalarIntNonNeg(value,field);          
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('    nupdate: [ Number of L-BFGS Hessian Updates to Store: {5} ] \n');
fprintf('      pgtol: [ Projected Gradient Tolerance: {1e-5} ] \n');
fprintf('\n');
