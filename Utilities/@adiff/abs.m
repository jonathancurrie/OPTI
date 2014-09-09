function c = abs(ad)
% ABS for adiff objects

if all(isreal(ad.x))&all(isreal(ad.dx(:)))
   c = adiff(abs(ad.x),rowmult(sign(ad.x), ad.dx),ad.root);
else
   c = sqrt(ad.*conj(ad));
end
