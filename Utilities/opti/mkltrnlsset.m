function options = mkltrnlsset(varargin)
%MKLTRNLSSET  Create or alter the options for Optimization with MKLTRNLS
%
% options = mkltrnlsset('param1',value1,'param2',value2,...) creates an MKLTRNLS
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the MKLTRNLS
% default.
%
% options = mkltrnlsset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = mkltrnlsset() creates an options structure with all fields set to
% mkltrnlsset defaults.
%
% mkltrnlsset() prints a list of all possible fields and their function.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'tolTR','tolSNGLR','tolSTEP','tolSTEPprec','trialStepSize'};
Defaults = {1e-6,1e-6,1e-6,1e-10,0.1}';        

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
    case {'toltr','tolsnglr','tolstep','tolstepprec'}
        err = opticheckval.checkScalarNonNeg(value,field);  
    case {'trialstepsize'}
        err = opticheckval.checkScalarBoundLLE(value,field,0,100);  
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('         tolTR: [ Trust Region Tolerance: {1e-6} ] \n');
fprintf('      tolSNGLR: [ Singular Convergence Tolerance: {1e-6} ] \n');
fprintf('       tolSTEP: [ Trial Step Tolerance: {1e-6} ] \n');
fprintf('   tolSTEPprec: [ Trial Step Precision: {1e-10} ] \n');
fprintf(' trialStepSize: [ Initial size of the trust region (boundary of trial step, 0.1 -> 100) {0.1} ] \n');
fprintf('\n');
