function ad1 = reciprocal(ad)
% RECIPROCAL computes the reciprocal of an adiff object,
% used for division.

ad1 = adiff(1./ad.x,rowmult(-1./ad.x.^2,ad.dx),ad.root);