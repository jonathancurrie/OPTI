function r = repmat(a,sz)
% REPMAT for adiff objects. sz must have all dimensions except
% the first equal to 1

if any(sz(2:end)~=1), error('repmat for adiff objects must have higher dimensions =1'); end
r = adiff( repmat(a.x,sz(1:2)), repmat(a.dx,sz(1:2)), a.root);