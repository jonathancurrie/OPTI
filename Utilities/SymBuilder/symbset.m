function options = symbset(varargin)
%SYMBSET  Create or alter the options for Optimization with SymBuilder
%
% options = symbset('param1',value1,'param2',value2,...) creates an
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the SymBuilder
% default.
%
% options = symbset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = symbset() creates an options structure with all fields set to
% symbset defaults.
%
% symbset() prints a list of all possible fields and their function.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

% Print out possible values of properties.
numberargs = nargin;
if (nargin == 0) && (nargout == 0)
    printfields();
    return
end
%Names and Defaults
Names = {'solver','display','cbmode','srckeep','use1stDerivs','use2ndDerivs','preallocate','optiOpts'}';
Defaults = {'auto','off','auto','no','yes','yes','yes',[]}';        

%Collect Sizes and lowercase matches         
m = size(Names,1); 
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

%If not using 1st derivatives, turn off 2nd too
if(strcmpi(options.use1stDerivs,'no')), options.use2ndDerivs = 'no'; end


function checkfield(field,value)
%Check a field contains correct data type

if isempty(value)
    return % empty matrix is always valid
end

switch lower(field)
    %String
    case {'use1stderivs','use2ndderivs','preallocate','srckeep'}
        if(ischar(value) && (strcmpi(value,'yes') || strcmpi(value,'no')))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a string with ''yes'' or ''no''',field);
        end
        
    case 'solver'
        if(ischar(value) && optiSolver(value))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a valid solver',field);
        end
    case 'cbmode'
        if(ischar(value) && any(strcmpi(value,{'auto','cppad','ccode','mcode'})))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be ''auto'', ''cppad'', ''ccode'', or ''mcode''',field);
        end
    case 'display'
        if(ischar(value) && any(strcmpi(value,{'off','iter','final'})))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be ''off'', ''iter'' or ''final''',field);
        end
    case {'solveropts','optiopts'}
        if(isstruct(value))
            valid = true;
        else
            valid = false; errmsg = sprintf('Parameter %s should be a structure',field);
        end
        
    otherwise  
        valid = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if valid 
    return;
else %error
    ME = MException('symbset:FieldError',errmsg);
    throwAsCaller(ME);
end


function printfields()
%Print out fields with defaults
fprintf('        solver: [ Optimization Solver {''auto''} ] \n');
fprintf('       display: [ Solver Display Level {''off''}, ''iter'', ''final'' ]\n');
fprintf('      optiOpts: [ OPTI options (via optiset) - overrides above {[]} ] \n');
fprintf('        cbmode: [ NL Problems: Callback Generation Mode {''auto''}, ''cppad'', ''ccode'', ''mcode'' ]\n');
fprintf('       srckeep: [ NL Problems: Keep C/C++ source generated by SymBuilder {''no''} ]\n');
fprintf('  use1stDerivs: [ NL Problems: Use first derivatives {''yes''} (Turns off use2ndDerivs too if no) ] \n');
fprintf('  use2ndDerivs: [ NL Problems: Use second derivatives {''yes''} ] \n');
fprintf('   preallocate: [ NL Problems: Preallocate return matrices/vectors {''yes''} ] \n');
fprintf('\n');
