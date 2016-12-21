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

%   Copyright (C) 2013 Jonathan Currie (IPL)

% Print out possible values of properties.
if ((nargin == 0) && (nargout == 0))
    printfields();
    return
end
%Names and Defaults
Names = {'scipopts','cipfile','gamsfile','testmode'}';
Defaults = {[],[],[],0}';        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar 0/1
    case 'testmode'
        err = opticheckval.checkScalar01(value,field);
    %char array
    case {'gamsfile','cipfile'}
        err = opticheckval.checkChar(value,field);    
    %cell array
    case 'scipopts'
        err = opticheckval.checkCell2Col(value,field); 
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('         scipopts: [ Set SCIP options using a 2D cell format: {''name1'', val1; ''name2'', val2; ...} (See Below) {[]} ] \n');...
fprintf('          cipfile: [ Write SCIP model to CIP file (will skip solving): {[]}, ''filename'' (Ensure Display is Off) ] \n');
fprintf('         gamsfile: [ Write SCIP model to GAMS file (will skip solving): {[]}, ''filename'' (Ensure Display is Off) ] \n');
fprintf('         testmode: [ Validate the nonlinear function generation (will skip solving): {0}, 1 ] \n');
fprintf('\n');

fprintf(' To set SCIP options not available via optiset, you may now set any available SCIP option using the ''scipopts'' field:\n');
fprintf(' - For a list of all available options, see http://scip.zib.de/doc/html/PARAMETERS.php\n');
fprintf(' - For example, to limit SCIP to the first feasible solution, use "scipset(''scipopts'',{''limits/solutions'',1})"\n');
fprintf(' - Note:\n    - Boolean, Int, LongInt and Real parameters are all assumed MATLAB doubles, and converted internally.\n');
fprintf('    - scipopts options take priority over all optiset options\n');
fprintf('\n');


