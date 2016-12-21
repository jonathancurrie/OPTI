function options = nloptset(varargin)
%NLOPTSET  Create or alter the options for Optimization with NLOPT
%
% options = nloptset('param1',value1,'param2',value2,...) creates an NLOPT
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the NLOPT
% default.
%
% options = nloptset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = nloptset() creates an options structure with all fields set to
% NLOPTSET defaults.
%
% nloptset() prints a list of all possible fields and their function.

%   Copyright (C) 2011 Jonathan Currie (IPL)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end

%Names and Defaults
Names = {'algorithm';'subalgorithm';'subtolrfun';'subtolafun';'submaxfeval';'submaxtime'};
Defaults = {'LN_AUGLAG';'LN_PRAXIS';1e-6;1e-6;1.5e3;1e3};         

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)   
    %Non-negative scalar double
    case {'subtolrfun','subtolafun'}
        err = opticheckval.checkScalarNonNeg(value,field); 
    %Non-zero scalar double
    case 'submaxtime'
        err = opticheckval.checkScalarGrtZ(value,field); 
    %Non-zero scalar integer
    case 'submaxfeval'  
        err = opticheckval.checkScalarIntGrtZ(value,field);  
    %NLOPT Algorithm
    case {'algorithm','subalgorithm'}
        err = opticheckval.checkNLOPTAlg(value,field);        
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults

fprintf('         algorithm: [ NLOPT Algorithm String ] \n');
fprintf('      subalgorithm: [ Suboptimizer NLOPT Algorithm String ] \n');
fprintf('        subtolrfun: [ Suboptimizer Relative Function Tolerance {1e-6} ] \n');
fprintf('        subtolafun: [ Suboptimizer Absolute Function Tolerance {1e-6} ] \n');
fprintf('       submaxfeval: [ Suboptimizer Maximum Function Evaluations {1.5e3} ] \n');
fprintf('        submaxtime: [ Suboptimizer Maximum Evaluation Time {1e3s} ] \n');
fprintf('\n');
