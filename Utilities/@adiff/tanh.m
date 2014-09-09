function c = tanh(a)
% TANH for adiff objects. 

c = adiff( tanh(a.x), rowmult(1./cosh(a.x).^2, a.dx), a.root);
