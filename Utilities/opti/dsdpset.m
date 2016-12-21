function options = dsdpset(varargin)
%DSDPSET  Create or alter the options for Optimization with DSDP
%
% options = dsdpset('param1',value1,'param2',value2,...) creates an DSDP
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the DSDP
% default.
%
% options = dsdpset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = dsdpset() creates an options structure with all fields set to
% dsdpset defaults.
%
% dsdpset() prints a list of all possible fields and their function.
%
% See supplied DSDP Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (IPL)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end
%Names and Defaults
Names = {'drho';'rpos';'r0';'penalty';'rho';'dbound';'gaptol';'rtol';'mu0';'maxtrust';...
         'steptol';'ptol';'pnormtol';'reuse';'zbar';'dlbound';'ybound';'fixed'};
Defaults = {1;0;-1;1e8;5;1e20;1e-7;1e-6;-1;1e10;0.05;1e-4;1e30;4;[];[];[];[]};        

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
    case {'r0','penalty','rho','mu0','dbound','pnormtol','zbar','dlbound','ybound'}
        err = opticheckval.checkScalarDbl(value,field);
    %Scalar non negative double
    case {'rtol','ptol','maxtrust','gaptol','steptol'}
        err = opticheckval.checkScalarNonNeg(value,field);      
    %Scalar 0/1
    case {'rpos','drho'}
        err = opticheckval.checkScalar01(value,field);   
    %Scalar non-negative integer
    case 'reuse'
        err = opticheckval.checkScalarIntNonNeg(value,field);
    %Double column vector
    case 'fixed'
        err = opticheckval.checkDualCol(value,field);        
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('          rpos: [ Use Penalty Parameter to enforce Feasibility {0}, 1 ] \n');
fprintf('          drho: [ Use Dynamic Strategy to choose Rho: 0, {1} ] \n');
fprintf('            r0: [ Set initial value for variable r in Dual {-1} ] \n');
fprintf('       penalty: [ Penalty Parameter Gamma {1e8} \n');
fprintf('           rho: [ Potential Parmater: {5} ] \n');
fprintf('           mu0: [ Barrier Parameter: {-1} ] \n');
fprintf('          rtol: [ Classify Dual as feasible if r is less than this tolerance: {1e-6} ] \n');
fprintf('          ptol: [ Classify Primal as feasible if infeasibility is less than this tolerance: {1e-4} ] \n');
fprintf('      maxtrust: [ Maximum Trust Radius on Step Direction: {1e10} ] \n');
fprintf('        gaptol: [ Convergence Gap Tolerance: {1e-7} ] \n');
fprintf('       steptol: [ Terminate solver if step length in Dual is below this tolerance: {0.05} ] \n');
fprintf('        dbound: [ Terminate solver if Dual Objective greater than this Value {1e20} ] \n');
fprintf('      pnormtol: [ Terminate the solver when duality gap is sufficiently small and pnorm is less than this quantity: {1e30} ] \n');
fprintf('         reuse: [ Number of times the Hessian of the Barrier Function will be reused {4} ] \n');
fprintf('          zbar: [ Upper bound on objective value at the solution {[]} ] \n');
fprintf('       dlbound: [ Lower bound on dual value {[]} ] \n');
fprintf('        ybound: [ Bound on the dual variables y {[]} ] \n');
fprintf('         fixed: [ Matrix of fixed variables, column 1 variable indices, column 2, fixed values {[]} ] \n');
fprintf('\n');
