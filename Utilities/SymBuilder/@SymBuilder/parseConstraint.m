function [str,cl,cu] = parseConstraint(str)
%PARSECON  Parse Constraint string and determine constraint type

%Replace double =
str = regexprep(str,'==','=');
%See if we have =
ind = strfind(str,'=');
if(isempty(ind))
    cl = 0;
    cu = 0; %assume = 0
else
    il = strfind(str,'<=');
    ig = strfind(str,'>=');
    try
        if(~isempty(il))
            cu = str2double(str(il+2:end));
            if(isnan(cu)) %assume variable
                cu = 0;
                str = [str(1:il-1) '-(' str(il+2:end) ')'];                                      
            else
                str = str(1:il-1);
                %Check if we can convert to a number (for ex 1<=1 we don't want)
                if(~isnan(str2double(str)))
                    str = []; cl = []; cu = [];
                    return;
                end
            end
            cl = -Inf;
        elseif(~isempty(ig))
            cl = str2double(str(ig+2:end));
            if(isnan(cl)) %assume variable
                cl = 0;
                str = [str(1:ig-1) '-(' str(ig+2:end) ')']; 
            else
                str = str(1:ig-1);
                %Check if we can convert to a number (for ex 1>=1 we don't want)
                if(~isnan(str2double(str)))
                    str = []; cl = []; cu = [];
                    return;
                end
            end
            cu = Inf;
        else 
            cl = str2double(str(ind+1:end));
            if(isnan(cl)) %assume variable
                cl = 0;
                str = [str(1:ind-1) '-(' str(ind+1:end) ')'];          
            else
                str = str(1:ind-1);
                %Check if we can convert to a number (for ex 1=1 we don't want)
                if(~isnan(str2double(str)))
                    str = []; cl = []; cu = [];
                    return;
                end
            end
            cu = cl;
        end
    catch ME
        error('There was an error parsing the RHS of the following equation:\n%s\n\nError: %s',str,ME.message);
    end
end