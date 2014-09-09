function c = exp(a)
% EXP for adiff objects. 

c = adiff(exp(a.x),rowmult(exp(a.x), a.dx),a.root);