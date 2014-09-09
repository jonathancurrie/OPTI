function c = ctranspose(ad)
% CTRANSPOSE returns the transpose of an adiff object.
% Unlike vectors, the transpose must be used with care, and
% is really only safe in places like bilinear or quadratic forms, 
% e.g. a'*B*c

c = adiff(ad.x',ad.dx,ad.root);
