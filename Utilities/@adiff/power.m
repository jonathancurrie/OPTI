function c = power(a,b)
% POWER implements a.^b, where either a or bor both is an adiff object.

switch [class(a),class(b)]
   
case 'adiffdouble'
   c = adiff( a.x.^b, rowmult(b.*a.x.^(b-1), a.dx), a.root);
   
case 'doubleadiff'
   c = adiff( a.^b.x, rowmult(a.^b.x.*log(a), b.dx), b.root);
   
case 'adiffadiff'
   c = exp(log(a).*b);
   % Commented out code does it the hard way:
   % cx  = a.x.^b.x;
   % cdx = rowmult(cx.*log(a.x), b.dx) + rowmult(a.x.^(b.x-1), a.dx);
   % c   = class(struct('x',cx,'dx',cdx),'adiff');
   
otherwise
   error(['Can''t do ',class(a),'.^',class(b)]);
   
end
