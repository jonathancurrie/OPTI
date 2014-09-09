function [var,lb,ub] = parseBound(str)
%PARSEBOUND  Parse bound and extract variable and bounds

%Strings to support
% lb <= m <= ub
% m >= lb
% m <= ub


%Find <=
ind = strfind(str,'<=');
if(isempty(ind))
    ind = strfind(str,'>=');
    if(isempty(ind))
        error('A bound must be of the form lb <= var <= ub, var >= lb, or var <= ub');
    else %lower bound
        var = strtrim(str(1:ind-1));
        lb = strtrim(str(ind+2:end));
        ub = 'Inf';
    end
else
    if(length(ind) > 1) %double bounds
        lb = strtrim(str(1:ind(1)-1));
        ub = strtrim(str(ind(2)+2:end));
        var = strtrim(str(ind(1)+2:ind(2)-1));
    else %upper bound
        var = strtrim(str(1:ind-1));
        ub = strtrim(str(ind+2:end));
        lb = '-Inf';
    end
end