function c = rdivide(a,b)
% RDIVIDE implements a./b, where either a or b is an adiff object.

switch [class(a),class(b)]
   
case 'adiffdouble'
   c = times(a,1./b);
   
case {'doubleadiff','adiffadiff'}
   c = times(a,reciprocal(b));
   
otherwise
   error(['Can''t divide (./) ',class(a),' and ',class(b)]);
   
end
