function options = mosekset(varargin)
%MOSEKSET  Create or alter the options for Optimization with MOSEK
%
% options = mosekset('param1',value1,'param2',value2,...) creates a MOSEK
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the MOSEK
% default.
%
% options = mosekset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = mosekset() creates an options structure with all fields set to
% mosekset defaults.
%
% mosekset() prints a list of all possible fields and their function.

%   Copyright (C) 2012 Jonathan Currie (IPL)

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end

%Names and Defaults
Names = {'maxiter';'ptol';'dtol';'maxbranch';'maxrelax';'tolabsgap';'tolrelgap';'tolint';'maxtime';'mskoption';'warnings';'display'};
Defaults = {400;1e-8;1e-8;-1;-1;0;1e-4;1e-5;-1;[];'critical';'off'};         

%Collect Sizes and lowercase matches         
m = size(Names,1); numberargs = nargin;
%Create structure with all names and default values
st = [Names,Defaults]'; options = struct(st{:});

% Check we have char or structure input. If structure, insert arguments
i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)
        break;
    end
    if ~isempty(arg)
        if ~isa(arg,'struct')
            error('An argument was supplied that wasn''t a string or struct!');
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

%Check we have even number of args left
if rem(numberargs-i+1,2) ~= 0
    error('You do not have a value specified for each field!');
end

%Go through each argument pair and assign to correct field
expectval = 0; %first arg is a name
while i <= numberargs
    arg = varargin{i};

    switch(expectval)
        case 0 %field
            if ~ischar(arg)
                error('Expected field name to be a string! - Argument %d',i);
            end
            j = find(strcmp(arg,Names) == 1);
            if isempty(j)  % if no matches
                error('Unrecognised parameter %s',arg);
            elseif(length(j) > 1)
                error('Ambiguous parameter %s',arg);
            end
            expectval = 1; %next arg is a value
        case 1
            checkfield(Names{j,:},lower(arg));
            options.(Names{j,:}) = lower(arg);
            expectval = 0;
    end
    i = i + 1;
end

if expectval %fallen off end somehow
    error('Missing value for %s',arg);
end



function checkfield(field,value)
%Check a field contains correct data type

if isempty(value)
    return % empty matrix is always valid
end

switch lower(field)    
    case {'maxiter','ptol','dtol','maxbranch','maxrelax','tolabsgap','tolrelgap','tolint','maxtime'} %general scalar
        if(isscalar(value) && isnumeric(value))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a scalar double',field);
        end      

    case 'mskoption' %structure
        if(isstruct(value))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a structure',field);
        end
        
    case 'warnings' %string
        if(ischar(value))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a char array',field);
        end                
    case 'display' %string
        switch(lower(value))
            case {'off','iter','final'}
                valid = true;
            otherwise
                valid = false;
                errmsg = sprintf('Parameter %s should be ''off'', ''iter'' or ''final'' ',field);
        end    
        
    otherwise  
        valid = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if valid 
    return;
else %error
    ME = MException('mosekset:FieldError',errmsg);
    throwAsCaller(ME);
end


function printfields()
%Print out fields with defaults

fprintf('  RELAXED SOLVER SETTINGS:\n');
fprintf('           maxiter: [ Maximum number of iterations {400} ] \n');       %MSK_IPAR_INTPNT_MAX_ITERATIONS
fprintf('              ptol: [ Primal Feasibility Tolerance {1e-8} ] \n');      %MSK_DPAR_INTPNT_TOL_PFEAS
fprintf('              dtol: [ Dual Feasibility Tolerance {1e-8} ] \n');        %MSK_DPAR_INTPNT_TOL_DFEAS

fprintf('  MIXED INTEGER SOLVER SETTINGS:\n');
fprintf('         maxbranch: [ Maximum number of branches {-1} ] \n');          %MSK_IPAR_MIO_MAX_NUM_BRANCHES
fprintf('          maxrelax: [ Maximum number of relaxations {-1} ] \n');       %MSK_IPAR_MIO_MAX_NUM_RELAXS
fprintf('         tolabsgap: [ Absolute optimality tolerance {0.0} ] \n');      %MSK_DPAR_MIO_TOL_ABS_GAP
fprintf('         tolrelgap: [ Relative optimality tolerance {1e-4} ] \n');     %MSK_DPAR_MIO_TOL_REL_GAP
fprintf('            tolint: [ Integer constraint tolerance {1e-5} ] \n');      %MSK_DPAR_MIO_TOL_ABS_RELAX_INT

fprintf('  GENERAL SETTINGS:\n');
fprintf('           maxtime: [ Maximum solve time {-1} ] \n');                  %MSK_DPAR_OPTIMIZER_MAX_TIME
fprintf('         mskoption: [ Structure of other MOSEK options {[]} ] \n');
fprintf('          warnings: [ {''on''} or ''off'' ] \n');
fprintf('           display: [ {''off''}, ''iter'', ''final'' ] \n');

fprintf('\nNote to use a MOSEK option other than in this list, pass it via mskoption, e.g.:\n');
fprintf('   >> mopt.MSK_DPAR_PRESOLVE_TOL_LIN_DEP = 1e-5;\n');
fprintf('   >> mosekset(''mskoption'',mopt);\n');
