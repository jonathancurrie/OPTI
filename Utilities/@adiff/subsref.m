function adsub = subsref(ad,s)
% SUBSREF picks out a part of an adiff object.
% The only available subsref is x(index) where x is an adiff object and index is
% some legal way of indexing a matlab vector

switch s(1).type
case '()' % subscripting of ad picks out just a part of it
   adsub = adiff( ad.x(s.subs{1}), ad.dx(s.subs{1},:), ad.root);
otherwise
   error('This kind of subsref isn''t allowed for adiff objects')
end

