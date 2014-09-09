function c = sinh(a)
% SINH for adiff objects. 

c = adiff( sinh(a.x), rowmult(cosh(a.x), a.dx), a.root);