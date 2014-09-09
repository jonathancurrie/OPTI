function c = times(a,b)
% TIMES implements a.*b, where either a or b is an adiff object.

switch [class(a),class(b)]
   
case 'adiffdouble'
   c = adiff(b.*a.x, rowmult(b, a.dx), a.root);
   
case 'doubleadiff'
   c = adiff(a.*b.x, rowmult(a, b.dx), b.root);
   
case 'adiffadiff'
   checkroot(a,b);
   c = adiff( a.x.*b.x, rowmult(b.x, a.dx)+rowmult(a.x, b.dx), a.root);
   
otherwise
   error(['Can''t multiply (.*) ',class(a),' and ',class(b)]);
   
end
