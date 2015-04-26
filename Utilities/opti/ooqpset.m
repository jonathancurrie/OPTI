function options = ooqpset(varargin)
%OOQPSET  Create or alter the options for Optimization with OOQP
%
% options = ooqpset('param1',value1,'param2',value2,...) creates an OOQP
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the OOQP
% default.
%
% options = ooqpset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = ooqpset() creates an options structure with all fields set to
% ooqpset defaults.
%
% ooqpset() prints a list of all possible fields and their function.
%
% See supplied OOQP Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'algorithm','linear_solver'};
Defaults = {'Gondzio','MA57'}';        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Misc String methods
    case 'algorithm'
        err = opticheckval.checkValidString(value,field,{'gondzio','mehrotra'});
    case 'linear_solver'
        err = opticheckval.checkValidString(value,field,{'pardiso','ma57','ma27'});
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('     algorithm: [ OOQP Solver Algorithm: {''Gondzio''}, ''Mehrotra'' ] \n');
fprintf(' linear_solver: [ OOQP Sparse Linear Solver: ''PARDISO'', {''MA57''}, ''MA27'' ] \n');
fprintf('\n');
