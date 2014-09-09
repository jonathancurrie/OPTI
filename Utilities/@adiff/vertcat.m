function c = vertcat(varargin)
% VERTCAT does vertical concatenation [a;b;...] of adiff objects and/or
% column vectors.
% The horizontal concatenation [a,b,...] is NOT supported

% Find the final dimensions of the result
nr = 0; nc = 0; issp = logical(0); root=[];
for i=1:length(varargin)
   if isa(varargin{i},'adiff')
      nr = nr+length(varargin{i});
      nvar = size(varargin{i}.dx,2);
      issp = issp&issparse(varargin{i});
      if isempty(root)|isempty(root.root), root = varargin{i}; end
      if nc==0, nc = nvar; end
      if nc~=nvar, error('All adiff objects must have the same number of variables'); end
   elseif isa(varargin{i},'double')
      if size(varargin{i},2)>1, error('Can only vertivally concatenate column vectors with adiff objects'); end
      nr = nr +length(varargin{i});
   end
end

% allocate the result
x = zeros(nr,1);
if issp
   dx=sparse(nr,nc);
else
   dx=zeros(nr,nc);
end

% fill the result
j = 1;
for i=1:length(varargin)
   nr = length(varargin{i});
   if isa(varargin{i},'adiff')
      checkroot(varargin{i},root);
      x(j:j+nr-1)=varargin{i}.x;
      dx(j:j+nr-1,:)=varargin{i}.dx;
   else
      x(j:j+nr-1)=varargin{i};
   end
   j = j+nr;
end

c = adiff( x, dx, root.root);


