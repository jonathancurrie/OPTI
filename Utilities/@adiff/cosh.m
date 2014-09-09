function c = cosh(a)
% COSH of an adiff object. 

c = adiff(cosh(a.x),rowmult(sinh(a.x), a.dx),a.root);
