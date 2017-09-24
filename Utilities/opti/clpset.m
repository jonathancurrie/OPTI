function options = clpset(varargin)
%CLPSET  Create or alter the options for Optimization with CLP
%
% options = clpset('param1',value1,'param2',value2,...) creates an CLP
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the CLP
% default.
%
% options = clpset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = clpset() creates an options structure with all fields set to
% clpset defaults.
%
% clpset() prints a list of all possible fields and their function.
%
% See supplied CLP Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (IPL)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'algorithm','primalTol','dualTol','doPresolve','numPresolvePasses','factorFreq','numberRefinements','primalObjLim','dualObjLim'};
Defaults = {'automatic',1e-7,1e-7,1,5,[],0,Inf,Inf};        

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
    case {'primalobjlim','dualobjlim'}
        err = opticheckval.checkScalarDbl(value,field);    
    %Scalar non zero double
    case {'primaltol','dualtol'}
        err = opticheckval.checkScalarGrtZ(value,field);    
    %Scalar non negative integer
    case {'numpresolvepasses','factorfreq','numberrefinements'}
        err = opticheckval.checkScalarIntNonNeg(value,field);          
    %Scalar integer with bounds
    case {'doscaling','dopresolve'}
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,1);    
    %Misc String methods
    case 'algorithm'
        err = opticheckval.checkValidString(value,field,{'DualSimplex','PrimalSimplex','PrimalSimplexOrSprint','Barrier','BarrierNoCross','Automatic'});
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('         algorithm: [ CLP Solver Algorithm: ''DualSimplex'', ''PrimalSimplex'', ''PrimalSimplexOrSprint'', ''Barrier'', ''BarrierNoCross'', {''Automatic''} ] \n');
fprintf('         primalTol: [ Primal Tolerance: {1e-7} ] \n');
fprintf('           dualTol: [ Dual Tolerance: {1e-7} ] \n');
fprintf('        doPresolve: [ Use CLP Presolver to Attempt to reduce Problem Size: Off (0), On {1} ] \n');
fprintf(' numPresolvePasses: [ Number of times to run the Presolver over the problem: {5} ] \n');
fprintf('        factorFreq: [ Simplex Factorization Frequency: {[]} (empty uses internal heuristic) ] \n');
fprintf(' numberRefinements: [ Number of iterative Simplex refinements: {0} ] \n');
fprintf('      primalObjLim: [ Primal Objective Limit: {Inf} ] \n');
fprintf('        dualObjLim: [ Dual Objective Limit: {Inf} ] \n');
% fprintf('        numThreads: [ Number of Cilk Worker Threads (LP Only and > 1 only with Aboca CLP Build): {1} ] \n');
% fprintf('          abcState: [ Aboca Partition Size: {[]} (empty uses internal heuristic) ] \n');
fprintf('\n');
