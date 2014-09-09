function c = mrdivide(a,b)
% MRDIVIDE implements a/b, where either a or b is an adiff object.
% It is mapped to a./b where possible (i.e. when b has length 1), and
% yields an error otherwise.

if prod(size(b))==1
   c = times(a,1./b);
else
   error('Cannot use / for adiff objects when denominator is nonscalar');
end
