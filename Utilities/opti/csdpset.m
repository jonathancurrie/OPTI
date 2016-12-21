function options = csdpset(varargin)
%CSDPSET  Create or alter the options for Optimization with CSDP
%
% options = csdpset('param1',value1,'param2',value2,...) creates an CSDP
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the CSDP
% default.
%
% options = csdpset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = csdpset() creates an options structure with all fields set to
% csdpset defaults.
%
% csdpset() prints a list of all possible fields and their function.
%
% See supplied CSDP Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (IPL)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'axtol';'atytol';'objtol';'pinftol';'dinftol';'minstepfrac';'maxstepfrac';...
         'minstepp';'minstepd';'usexzgap';'tweakgap';'affine';'perturbobj';'fastmode';...
         'objconstant';'writeprob';'writesol'};
Defaults = {1e-8,1e-8,1e-8,1e8,1e8,0.9,0.97,1e-8,1e-8,1,0,0,1,0,0,[],[]}';        

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
    case {'perturbobj','objconstant'}
        err = opticheckval.checkScalarDbl(value,field);
    %Scalar non negative double
    case {'axtol','atytol','objtol','pinftol','dinftol','minstepfrac','maxstepfrac','minstepp','minstepd'}
        err = opticheckval.checkScalarNonNeg(value,field);      
    %Scalar 0/1
    case {'usexzgap','tweakgap','affine','fastmode'}
        err = opticheckval.checkScalar01(value,field);   
    %Char Array
    case {'writeprob','writesol'}
        err = opticheckval.checkChar(value,field);           
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('         axtol: [ Primal Feasibility Tolerance: {1e-8} ] \n');
fprintf('        atytol: [ Dual Feasibility Tolerance: {1e-8} ] \n');
fprintf('        objtol: [ Relative Duality Gap Tolerance: {1e-8} ] \n');
fprintf('       pinftol: [ Primal Infeasibility Tolerance: {1e8} \n');
fprintf('       dinftol: [ Dual Infeasibility Tolerance: {1e8} ] \n');
fprintf('   minstepfrac: [ Minimum step fraction to edge of feasible region: {0.9} ] \n');
fprintf('   maxstepfrac: [ Maximum step fraction to edge of feasible region: {0.97} ] \n');
fprintf('      minstepp: [ Primal line search minimum step size before failure: {1e-8} ] \n');
fprintf('      minstepd: [ Dual line search minimum step size before failure: {1e-8} ] \n');
fprintf('      usexzgap: [ Use objective function duality gap instead of tr(XZ) gap: 0, {1} ] \n');
fprintf('      tweakgap: [ Attempt to fix negative duality gaps if usexzgap = 0: {0}, 1 ] \n');
fprintf('        affine: [ Only take primal-dual affine steps (do not use barrier term): {0}, 1 ] \n');
fprintf('    perturbobj: [ Level of objective function perturbation (useful in unbounded problems): {1.0} ] \n');
fprintf('      fastmode: [ Sacrifice accuracy for faster execution: {0}, 1 ] \n');
fprintf('   objconstant: [ Constant value to add to the objective: {0.0} ] \n');
fprintf('     writeprob: [ Filename (including path) to write the entered problem to in SDPA sparse format: {[]} ] \n');
fprintf('      writesol: [ Filename (including path) to write the solution to in SDPA sparse format: {[]} ] \n');
fprintf('\n');
