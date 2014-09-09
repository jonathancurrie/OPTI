function c = log10(a)
% LOG10 for adiff objects. 

c = adiff(log(a.x)/log(10), rowmult(1./a.x/log(10), a.dx), a.root);