function c = conj(a)
% CONJ for adiff objects
c=adiff(conj(a.x),conj(a.dx),a.root);