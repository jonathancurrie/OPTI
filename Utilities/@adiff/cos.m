function c = cos(a)
% COS of an adiff object. 

c = adiff(cos(a.x),rowmult(-sin(a.x), a.dx),a.root);
