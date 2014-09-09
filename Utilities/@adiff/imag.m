function c = imag(a)
% IMAG for adiff objects
c=adiff(imag(a.x),imag(a.dx),a.root);