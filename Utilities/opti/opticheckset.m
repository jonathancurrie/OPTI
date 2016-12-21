function options = opticheckset(Names,Defaults,checkfun,varargin)
%OPTICHECKSET  Default Algorithm for dealing with solver option setting

%   Copyright (C) 2013 Jonathan Currie (IPL)

%Transpose as required
if(size(Names,2) > 1), Names = Names'; end
if(size(Defaults,2) > 1), Defaults = Defaults'; end

%Collect Sizes and lowercase matches         
m = size(Names,1); 
numberargs = numel(varargin);
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
                checkfun(Names{j,:},val);
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
                error('Expected fieldname (argument %d) to be a string!',i);
            end
            j = find(strcmp(arg,Names) == 1);
            if isempty(j)  % if no matches
                error('Unrecognised parameter ''%s''',arg);
            elseif(length(j) > 1)
                error('Ambiguous parameter ''%s''',arg);
            end
            expectval = 1; %next arg is a value
        case 1
            if(~isempty(arg))   
                if(ischar(arg)), arg = lower(arg); end
                checkfun(Names{j,:},arg);
                options.(Names{j,:}) = arg;
            end
            expectval = 0;
    end
    i = i + 1;
end

if expectval %fallen off end somehow
    error('Missing value for final argument ''%s''',arg);
end