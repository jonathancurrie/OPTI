function c = uminus(ad)
% UMINUS implements unary minus for adiff objects

% CHECKED

c = adiff( -ad.x, -ad.dx, ad.root);