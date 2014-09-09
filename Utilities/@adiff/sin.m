function c = sin(a)
% SIN for adiff objects. 

c = adiff( sin(a.x), rowmult(cos(a.x), a.dx), a.root);