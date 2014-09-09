function options = pswarmset(varargin)
%PSWARMSET  Create or alter the options for Optimization with PSwarm
%
% options = pswarmset('param1',value1,'param2',value2,...) creates an
% PSWARM options structure with the parameters 'param' set to their 
% corresponding values in 'value'. Parameters not specified will be set to 
% the PSWARM default.
%
% options = pswarmset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = pswarmset() creates an options structure with all fields set to
% PSWARMSET defaults.
%
% pswarmset() prints a list of all possible fields and their function.

%   Copyright (C) 2012 Jonathan Currie (I2C2)

% Print out possible values of properties.
if ((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'swarm_size';'vectorized';'mu';'nu';'iweight';'fweight';'delta';'idelta';'ddelta'};
Defaults = {42;0;0.5;0.5;0.9;0.4;Inf;2.0;0.5};         

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar double
    case {'mu','nu','iweight','fweight','delta','idelta','ddelta'}
        err = opticheckval.checkScalarDbl(value,field);      
    %Scalar 0/1
    case 'vectorized'
        err = opticheckval.checkScalar01(value,field);   
    %Non-zero scalar double
    case 'swarm_size'
        err = opticheckval.checkScalarIntGrtZ(value,field);   
        
    otherwise  
        err = sprintf('Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('        swarm_size: [ Swarm Size {42} ] \n');
fprintf('        vectorized: [ Objective function is vectorized {0}, 1 ] \n');
fprintf('                mu: [ Cognitial Parameter {0.5} ] \n');
fprintf('                nu: [ Social Parameter {0.5} ] \n');
fprintf('           iweight: [ Initial Weight {0.9} ] \n');
fprintf('           fweight: [ Final Weight {0.4} ] \n');
fprintf('             delta: [ Initial Delta {Inf} ] \n');
fprintf('            idelta: [ Increase Delta {2.0} ] \n');
fprintf('            ddelta: [ Decrease Delta {0.5} ] \n');
fprintf('\n');

