function c = tan(a)
% TAN for adiff objects. 

c = adiff( tan(a.x), rowmult(1./cos(a.x).^2, a.dx), a.root);
