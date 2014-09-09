function fnames = nestedfieldnames(field,pre)
%Return a cell array of strings containing nested structure args(recursive)
fnames = {};
if(~isstruct(field))
    error('This function only works on structures');
end
fn = fieldnames(field);
no = length(fn);
if(no < 1)
    error('The supplied structure must have at least one field');
end
if(nargin < 2)
    pre = '';
else
    pre = [pre '.'];
end

for i = 1:no
    if(isstruct(field.(fn{i})))
        fnames{i} = nestedfieldnames(field.(fn{i}),[pre fn{i}]); %#ok<AGROW>
    else
        fnames{i} = [pre fn{i}]; %#ok<AGROW>
    end
end
%Flatten cell array
fnames = flatten(fnames);

function C = flatten(A)
% Ref: http://www.mathworks.com/matlabcentral/fileexchange/27009-flatten-nested-cell-arrays/content/flatten.m
% (recursive)
C = {};
for i=1:numel(A)  
    if(~iscell(A{i}))
        C = [C,A{i}]; %#ok<AGROW>
    else
        Ctemp = flatten(A{i});
        C = [C,Ctemp{:}]; %#ok<AGROW>
    end
end