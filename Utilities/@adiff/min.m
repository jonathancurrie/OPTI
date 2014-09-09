function c = min(a,b)
% MIN for adiff objects. This selects the value of a which is a maximum. 
% Both min(a) and min(a,b) will work. min(a,[],dim) is meaningless because
% the adiff object is a column vector; use min(a) instead.

if nargin==1
   [y,i] = min(a.x);
   c = adiff(y,a.dx(i,:),a.root);
else
   switch [class(a),class(b)]
   case 'doubleadiff'
      ok = b.x<a;
      c = adiff(min(a,b.x), rowmult(ok,b.dx), b.root);
   case 'adiffdouble'
      ok = a.x<b;
      c = adiff(min(a.x,b), rowmult(ok,a.dx), a.root);
   case 'adiffadiff'
      checkroot(a,b);
      if size(a.dx,1)~=size(b.dx,1)
         if size(a.dx,1)==1
            a.dx = repmat(a.dx,size(b.dx,1),1);
         elseif size(b.dx,1)==1
            b.dx = repmat(b.dx,size(a.dx,1),1);
         end
      end
      ok = a.x<b.x;
      c = adiff(min(a.x,b.x), rowmult(ok,a.dx)+rowmult(1-ok,b.dx), a.root);
   otherwise
      error(['Can''t compute max of ',class(a),' and ',class(b)]);
   end
end
