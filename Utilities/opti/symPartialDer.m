function symder = symPartialDer(symfun,var,ncol,ind,der2)
%SYMPARTIALDER  Solve for the Partial Derivative of a Supplied Symbolic Function
%
%   jac = symPartialDer(fun,var) uses the symbolic toolbox to 
%   automatically generate the partial derivative matrix of the symbolic
%   function fun. var specifies the variable to differentiate with respect 
%   to (as a string).

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 5), der2 = false; end
if(nargin < 4), ind = []; end
if(nargin < 3), ncol = 0; end

%Find unique indices
ind = unique(ind);
%Check we have enough vars, otherwise manually generate the jacobian/hessian
if(ncol && (length(ind) ~= ncol))
    symder = '';
    if(der2)
        %Manually Generate Hessian
        for i = 1:ncol
            for j = 1:ncol
                symder = [symder diff(diff(symfun,sym(sprintf('%s%d',var,i))),sym(sprintf('%s%d',var,j)))]; %#ok<AGROW>
            end
        end    
        symder = reshape(symder,ncol,ncol);
    else
        %Manually Generate Jacobian        
        for i = 1:ncol
            symder = [symder diff(symfun,sym(sprintf('%s%d',var,i)))]; %#ok<AGROW>
        end    
    end
else
    %Determine variables contained within symbolic function
    vars = sort(findSymVars(symfun,var));
    if(der2)
        %Calculate partial derivative Hessian
        symder = hessian(symfun,vars);
    else
        %Calculate partial derivative jacobian
        symder = jacobian(symfun,vars);
    end
end


function vars = findSymVars(symfun,var)
str = char(symvar(symfun));
%Remove 'matrix', if present
if(~isempty(strfind(str,'matrix')))
    str = str(7:end-1);
end
ind = strfind(str,var);
if(~isempty(ind))
    %If no brackets, then only one variable
    if(isempty(strfind(str,',')) && isempty(strfind(str,']')))
        vars = symvar(symfun);
    else
        %Symvar always puts in alpha order, use to our advantage, find next comma or bracket from last index
        for i = ind(end)+1:length(str)
            if(any(strcmp(str(i),{',',']'})))
                break;
            end
        end
        if(i == length(str))
            error('Didn''t find the end of the symbolic variable string??');
        end
        vars = sym(sprintf('[%s]',str(ind(1):i-1)));
    end
else
    vars = [];
end