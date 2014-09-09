function c = log(a)
% LOG for adiff objects. 

c = adiff(log(a.x),rowmult(1./a.x, a.dx),a.root);