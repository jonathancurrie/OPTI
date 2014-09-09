function options = scipset(varargin)
%SCIPSET  Create or alter the options for Optimization with SCIP
%
% options = scipset('param1',value1,'param2',value2,...) creates an SCIP
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the SCIP
% default.
%
% options = scipset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = scipset() creates an options structure with all fields set to
% SCIPSET defaults.
%
% scipset() prints a list of all possible fields and their function.
%
% See supplied SCIP Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

% Print out possible values of properties.
if ((nargin == 0) && (nargout == 0))
    printfields();
    return
end
%Names and Defaults
Names = {'gamsfile','cipfile','maxpresolve','testmode'}';
Defaults = {[],[],-1,0}';        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Integer -1 <= value < Inf
    case 'maxpresolve'
        err = opticheckval.checkScalarIntBoundLEL(value,field,-1,Inf);
    %Scalar 0/1
    case 'testmode'
        err = opticheckval.checkScalar01(value,field);
    %char array
    case {'gamsfile','cipfile'}
        err = opticheckval.checkChar(value,field);         
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('         gamsfile: [ Write SCIP model to GAMS file (will skip solving): {[]}, ''filename'' ] \n');
fprintf('         testmode: [ Validate the nonlinear function generation (will skip solving): {0}, 1 ] \n');
fprintf('\n');
