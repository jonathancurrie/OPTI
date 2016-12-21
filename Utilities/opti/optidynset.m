function options = optidynset(varargin)
%OPTIDYNSET  Create or alter the options for Dynamic Optimization with OPTI
%
% options = optidynset('param1',value1,'param2',value2,...) creates an OPTI
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the OPTI
% default.
%
% options = optidynset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = optidynset() creates an options structure with all fields set to
% OPTI defaults.
%
% optidynset() prints a list of all possible fields and their function.

%   Copyright (C) 2013 Jonathan Currie (IPL)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'integrator','sensitivity','dfdz','dfdp','stateIndex','initialT','odeMaxTime','odeOpts'};
Defaults = {'ode45',[],[],[],[],[],30,[]};         

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end

function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Function Handle
    case {'dfdz','dfdp'}
        err = opticheckval.checkFunHandle(value,field);
    %Logical Vector or Double Vector
    case 'stateindex'
        err = opticheckval.checkDblVecOrLogVec(value,field);
    %Double >= 0
    case {'initialt','odemaxtime'}
        err = opticheckval.checkScalarBoundLEL(value,field,0,inf);
    %Struct
    case 'odeopts'
        err = opticheckval.checkStruct(value,field);
    %Other Strings
    case 'integrator'
        err = opticheckval.checkValidString(value,field,{'ode45','ode15s','ode23','ode113','ode23t','ode23tb','ode23s','ode15i'});
    case 'sensitivity'        
        err = opticheckval.checkValidString(value,field,{'cs','nd','ad','user','none'});
        
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults

fprintf('\n OPTIDYNSET Fields:\n');
fprintf('       integrator: [ ODE Integrator: {''ode45''}, ''ode15s'', ''ode23'', ''ode113'', ''ode23t'', ''ode23tb'', ''ode23s'', ''ode15i'' ] \n');
fprintf('      sensitivity: [ ODE Sensitivity Derivatives (defaults to User if supplied, otherwise ND): Complex Step ''CS'', Numerical ''ND'', Automatic ''AD'', User Supplied ''User'', Don''t Use Sensitivity ''None'' ] \n');
fprintf('             dfdz: [ Jacobian of ODE with Respect to State Variables (z) ] \n');
fprintf('             dfdp: [ Jacobian of ODE with Respect to Parameters (p) ] \n');
fprintf('       stateIndex: [ Indices of States to Fit to Measured Data ] \n');
fprintf('         initialT: [ Initial Time to Start ODE Integrator {[]} ] \n');
fprintf('       odeMaxTime: [ Maximum Time the ODE Integrator Can Run per Call [sec] {30} ] \n');
fprintf('          odeOpts: [ ODE Integrator Specific Options ] \n');


