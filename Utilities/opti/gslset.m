function options = gslset(varargin)
%GSLSET  Create or alter the options for Optimization with GSL
%
% options = gslset('param1',value1,'param2',value2,...) creates a GSL
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the GSL
% default.
%
% options = gslset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = gslset() creates an options structure with all fields set to
% gslset defaults.
%
% gslset() prints a list of all possible fields and their function.

%   Copyright (C) 2017 Jonathan Currie (IPL)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'factorUp','factorDown','avmax','h_df','trustRegionSolver', 'scalingMethod', 'linearSolver'}';
Defaults = {3.0,2.0,0.75,sqrt(eps),'lm','more','qr'}';        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar double > 0
    case {'factorup','factordown','avmax','h_df'}
        err = opticheckval.checkScalarGrtZ(value,field);  
    case {'trustregionsolver'}
        err = opticheckval.checkValidString(value,field, {'lm','lmaccel','dogleg','ddogleg','subspace2D'});  
    case {'scalingmethod'}
        err = opticheckval.checkValidString(value,field, {'levenberg','marquardt','more'}); 
    case {'linearsolver'}
        err = opticheckval.checkValidString(value,field, {'qr','cholesky','svd'}); 
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults

fprintf('\n GSL Nonlinear Least Squares options:\n');
fprintf('          factorUp: [ Factor for Increasing Trust Radius: {3.0} ] \n');
fprintf('        factorDown: [ Factor for Decreasing Trust Radius: {2.0} ] \n');
fprintf('             avmax: [ Max allowed |a|/|v|: {0.75} ] \n');
fprintf('              h_df: [ Step-size for Finite Difference Jacobian (only if OPTI class is not used): {sqrt(eps)} ] \n');
fprintf(' trustRegionSolver: [ Trust Region Suproblem Method: {''lm''}, ''lmaccel'', ''dogleg'', ''ddogleg'', ''subspace2D'' ] \n');
fprintf('     scalingMethod: [ Scaling Method: {''more''}, ''levenberg'', ''marquardt'' ] \n');
fprintf('      linearSolver: [ Solver Method: {''qr''}, ''cholesky'', ''svd'' ] \n');
fprintf('\n');
