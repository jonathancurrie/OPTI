function [n,m] = size(ad,dim)
% SIZE for adiff objects

[n,m] = size(ad.x);
if nargin==2 & dim==2, n=m; end
if nargout<2, n = [n,m]; end

