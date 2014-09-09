function c = mtimes(a,b)
% MTIMES implements a*b, where one of a or b is an adiff object.
% If both are adiff objects, this will only work if at least one
% has unit length, or a is a transposed adiff object. Otherwise there
% will be an error of some sort.

if all(size(a)==1)|all(size(b)==1)
   c = times(a,b);
   return
end

switch [class(a),class(b)]
   
case 'adiffdouble'
   c = adiff(a.x*b,b'*a.dx, a.root);
    
case 'doubleadiff'
   c = adiff(a*b.x,a*b.dx,b.root);
   
case 'adiffadiff'
   % this will only work if a was transposed
   checkroot(a,b);
   c = adiff(a.x*b.x, a.x*b.dx+b.x'*a.dx, a.root);
   
otherwise
   error(['Can''t matrix multiply ',class(a),' and ',class(b)]);
   
end
