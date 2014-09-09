function c = mldivide(a,b)
% MLDIVIDE implements a\b, where a is double and b is an adiff object, or
% when a is an adiff object with length 1.

switch [class(a),class(b)]
case 'doubleadiff'
   c = adiff(a\b.x, a\b.dx, b.root);
case {'adiffdouble','adiffadiff'}
   if length(a)==1
      c = times(b,1./a);
   else
      error('Cannot do a\b when a is an adiff object with length>1')
   end
otherwise 
   error(['\ not defined for ',class(a),' and ',class(b)]);
end
