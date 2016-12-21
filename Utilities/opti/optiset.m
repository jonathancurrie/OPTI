function options = optiset(varargin)
%OPTISET  Create or alter the options for Optimization with OPTI
%
% options = optiset('param1',value1,'param2',value2,...) creates an OPTI
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the OPTI
% default.
%
% options = optiset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = optiset() creates an options structure with all fields set to
% OPTI defaults.
%
% optiset() prints a list of all possible fields and their function.

%   Copyright (C) 2011 Jonathan Currie (IPL)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'solver';'maxiter';'maxfeval';'maxnodes';'maxtime';'tolrfun';'tolafun';'tolint';'solverOpts';'dynamicOpts';'iterfun';'warnings';'display';'derivCheck'};
Defaults = {'auto';1500;1e4;1e4;1000;1e-7;1e-7;1e-5;[];[];[];'critical';'off';'off'};         

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end

%Ensure backwards compatible
if(strcmpi(options.warnings,'on'))
    options.warnings = 'all';
elseif(strcmpi(options.warnings,'off'))
    options.warnings = 'none';
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar non negative double
    case {'tolafun','tolrfun','tolint'}
        err = opticheckval.checkScalarNonNeg(value,field);    
    %Scalar non zero double
    case 'maxtime'
        err = opticheckval.checkScalarGrtZ(value,field);  
    %Scalar non zero integer
    case {'maxiter','maxfeval','maxnodes'}
        err = opticheckval.checkScalarIntGrtZ(value,field);
    %Struct
    case {'solveropts','dynamicopts'}
        err = opticheckval.checkStruct(value,field);
    %Function Handle
    case 'iterfun'
        err = opticheckval.checkFunHandle(value,field);
    %Other misc
    case 'display'
        err = opticheckval.checkValidString(value,field,{'off','iter','final'});
    case 'warnings'        
        err = opticheckval.checkValidString(value,field,{'on','critical','off','all','none'});
    case 'derivcheck'        
        err = opticheckval.checkValidString(value,field,{'on','off'});
    case 'solver'
        err = opticheckval.checkOptiSolver(value,field);
        
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults

solvers = [{'AUTO'} optiSolver('all')];
len = length(solvers);
str = '';
for i = 1:len
    if(i < len)
        str = [str solvers{i} ', ']; %#ok<AGROW>
    else
        str = [str solvers{i}]; %#ok<AGROW>
    end
end   
str = regexprep(str,'AUTO','{AUTO}');

fprintf('\n OPTISET Fields:\n');
fprintf(['            solver: [ ' str ' ] \n']);
fprintf('           maxiter: [ Maximum Solver Iterations {1.5e3} ] \n');
fprintf('          maxfeval: [ Maximum Function Evaluations {1e4} ] \n');
fprintf('          maxnodes: [ Maximum Integer Solver Nodes {1e4} ] \n');
fprintf('           maxtime: [ Maximum Solver Evaluation Time {1e3s} ] \n');
fprintf('           tolrfun: [ Relative Function Tolerance {1e-7} ] \n');
fprintf('           tolafun: [ Absolute Function Tolerance {1e-7} ] \n');
fprintf('            tolint: [ Absolute Integer Tolerance {1e-5} ] \n');
fprintf('        solverOpts: [ Solver Specific Options Structure ] \n');
fprintf('       dynamicOpts: [ Dynamic Optimization Options Structure ] \n');
fprintf('           iterfun: [ Iteration Callback Function, stop = iterfun(iter,fval,x) {} ] \n');
fprintf('          warnings: [ ''all'' or {''critical''} or ''none'' ] \n');
fprintf('           display: [ {''off''}, ''iter'', ''final'' ] \n');
fprintf('        derivCheck: [ Derivative Checker {''off''}, ''on'' ] \n');

