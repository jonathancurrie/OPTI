function symder = symPartialDer(symfun,var,ncol,ind,der2)
%SYMPARTIALDER  Solve for the Partial Derivative of a Supplied Symbolic Function
%
%   jac = symPartialDer(fun,var) uses the symbolic toolbox to 
%   automatically generate the partial derivative matrix of the symbolic
%   function fun. var specifies the variable to differentiate with respect 
%   to (as a string).

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

if(nargin < 5), der2 = false; end
if(nargin < 4), ind = []; end
if(nargin < 3), ncol = 0; end

haveStr2Sym = (exist('str2sym','file') == 2);

%Find unique indices
ind = unique(ind);
%Check we have enough vars, otherwise manually generate the jacobian/hessian
if(ncol && (length(ind) ~= ncol))
    if (haveStr2Sym)
        symder = str2sym('');
    else
        symder = '';
    end
    if(der2)
        %Manually Generate Hessian
        for i = 1:ncol
            var1str = sprintf('%s%d',var,i);
            if (haveStr2Sym)
                var1 = str2sym(var1str);
            else
                var1 = sym(var1str);
            end
            for j = 1:ncol
                var2str = sprintf('%s%d',var,j);
                if (haveStr2Sym)
                    var2 = str2sym(var2str);
                else
                    var2 = sym(var2str);
                end
                symder = [symder diff(diff(symfun,var1),var2)]; %#ok<AGROW>
            end
        end    
        symder = reshape(symder,ncol,ncol);
    else
        %Manually Generate Jacobian        
        for i = 1:ncol
            var1str = sprintf('%s%d',var,i);
            if (haveStr2Sym)
                var1 = str2sym(var1str);
            else
                var1 = sym(var1str);
            end
            symder = [symder diff(symfun,var1)]; %#ok<AGROW>
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
        if(i == length(str) && str(i) ~= ']')
            error('Didn''t find the end of the symbolic variable string??');
        end
        varstr = sprintf('[%s]',str(ind(1):i-1));
        if (exist('str2sym','file') == 2)
            vars = str2sym(varstr);
        else
            vars = sym(varstr);
        end
    end
else
    vars = [];
end