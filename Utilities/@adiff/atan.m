function c = atan(a)
% ATAN for adiff objects. 

c = adiff(atan(a.x),rowmult(1./(1+a.x.^2), a.dx),a.root);
