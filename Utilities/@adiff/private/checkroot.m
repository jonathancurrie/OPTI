function checkroot(a,b)
% CHECKROOT checks that the root members of a&b are the same

if isempty(a.root)|isempty(b.root), return, end
if (length(a.root)~=length(b.root))|any(a.root~=b.root)
   error('adiff objects must come from the same root')
end
