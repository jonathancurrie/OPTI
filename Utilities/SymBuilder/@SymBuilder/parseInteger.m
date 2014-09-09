function [var,xtype] = parseInteger(str)
%PARSEBOUND  Parse integer and extract variable and type

%Strings to support
% var = I
% var = B
% var = C


%Find =
ind = strfind(str,'=');
if(isempty(ind))
    error('An integer declaration must be of the form var = I/B/C');
else
    var = strtrim(str(1:ind-1));
    xx = upper(strtrim(str(ind+1:end)));
    switch(xx)
        case {'C','I','B'}
            xtype = xx;
        otherwise
            error('Unknown integer type: %s',xx);
    end
end