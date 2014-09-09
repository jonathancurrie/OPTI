function da = diff(ad,n)
% DIFF for adiff objects, returns the discrete derivative
% along the value of the adiff object.

if nargin<2, n=1; end
da = adiff(diff(ad.x,n),diff(ad.dx,n),ad.root);
