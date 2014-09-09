function a = subsasgn(a,s,b)
% SUBSASGN for adiff objects

if length(s)>1, error('Subscripted assignment too complicated'); end

if isempty(a)
   if issparse(b)
      a = adiff([], sparse(0,size(b.dx,2)), []);
   else
      a = adiff( [], zeros(0,size(b.dx,2)), []);
   end
end

if length(s.subs)==1, s.subs{end+1}=':'; end

switch [class(a),class(b)]
case 'doubleadiff'
   a = subsasgn(a,s,b.x);
case 'adiffdouble'
   a.x  = subsasgn(a.x,s,b);
   if isempty(b)
      a.dx = subsasgn(a.dx,s,[]);
   else
      a.dx = subsasgn(a.dx,s,0);
   end      
case 'adiffadiff'
   checkroot(a,b);
   a.x  = subsasgn(a.x,s,b.x);
   a.dx = subsasgn(full(a.dx),s,full(b.dx));
otherwise
   error(['cannot assign ',class(b),' to ',class(a)]);
end