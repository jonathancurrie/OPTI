function a = rowmult(v,m)
% ROWMULT(v,m) implements diag(v)*m where v is a column vector (efficiently?)

% This is a good candidate for a mex file.

if size(v,2)>1, error('rowmult: can only deal with column vectors'); end

% There are three separate cases to deal with
if size(m,1)==1
   a = v(:)*m;
elseif length(v)==1
   a = v*m;
elseif length(v)==size(m,1)
   if issparse(m)
      [r,c,mv] = find(m);
      [nr,nc]  = size(m);
      a = sparse(r,c,mv.*v(r),nr,nc);
   else 
      a  = repmat(v(:),[1,size(m,2)]).*m;
   end
else
   error('Bad sizes in rowmult')
end
