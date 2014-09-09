function c = asin(a)
% ASIN for adiff objects. 

c = adiff(asin(a.x),rowmult(1./sqrt(1-a.x.^2), a.dx),a.root);