function fun = sym2fun(sobj,svar,var)
% SYM2FUN  Convert Symbolic Expression into Matlab Function Handle    

if(nargin < 3 || isempty(var)), var = 'x'; end

%Build cell array to store indexed vars
ivar = convIndex(svar);

%Subs out individual symbolic variables into our indexed list
eq = subs(sobj,svar,ivar);

%Build a MATLAB function
str = char(eq);
%Remove 'matrix' if is found
ind = strfind(str,'matrix');
if(ind)
    str = str(ind+7:end-1);
end
%Ensure we have right size
[r,c] = size(eq);
if(c==1 && r > 1)
    str = regexprep(str,',',';'); %force column
elseif(c>1 && r>1)
    str = regexprep(str,'],','];'); %keep matrix
end  
%Build function handle
fun = str2func(['@(' var ') ' str]);


%Convert index from m1 to m(1)
function [ivar,var] = convIndex(svar)            
ivar = cell(size(svar));
%Find indexing number
for i = 1:length(svar)
    str = char(svar(i));
    for j = 1:length(str)
        if(~isletter(str(j)))
            break;
        end
    end
    if(i == 1) %assume first variable contains decision variable
        var = str(1:j-1);
    end
%   ivar{i} = sprintf('%s(%s)',str(1:j-1),str(j:end)); %this uses user index...
    ivar{i} = sprintf('%s(%d)',str(1:j-1),i); %this uses index based on array
end