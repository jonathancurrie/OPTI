function [group,name] = parseResExp(str)
%PARSERESEXP  Parse result expression and extract group and name

%Find =
ind = strfind(str,':');
if(isempty(ind))
    error('An result expression name must be of the form group:name');
else
    name = strtrim(str(ind+1:end));
    group = strtrim(str(1:ind-1));
end