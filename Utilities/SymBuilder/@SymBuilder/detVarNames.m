function varn = detVarNames(svar)
%Determine variable names, in order, index free
%e.g. [ h1, h2, m1, m2, m3] becomes {h [1 2]; m [3 4 5]}

%Initialize
varn = []; k = 1;
%Search each var
for i = 1:length(svar)
    str = char(svar(i)); found = 0;
    for j = 1:length(str)
        if(uint8(str(j)) < 65 || uint8(str(j)) > 126)
            found = 1;
            break; %found a number, assume end of name
        end
    end
    %If we found it, remove number (last index found)
    if(found)        
        nvar = str(1:j-1); %New variable name
    else
        nvar = str;
    end
    %See where to place it
    if(isempty(varn))
        varn = {nvar i};
        k = k + 1;
    elseif(strcmp(varn{k-1,1},nvar)) %same as old one
        varn{k-1,2} = [varn{k-1,2} i]; %#ok<AGROW>
    else %new one
        varn = [varn; {nvar} i]; %#ok<AGROW>
        k = k + 1;
    end
end