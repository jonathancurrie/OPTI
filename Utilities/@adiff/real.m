function c = real(a)
% REAL for adiff objects
c=adiff(real(a.x),real(a.dx),a.root);