function [str,exp] = parseExpression(str)
%PARSEEXPRESSION  Parse expression and extract variable and expression

%Find =
ind = strfind(str,'=');
if(isempty(ind))
    error('An expression must be of the form a = b*x');
else
    exp = strtrim(str(ind+1:end));
    str = strtrim(str(1:ind-1));
end