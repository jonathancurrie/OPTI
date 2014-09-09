function [x,dx] = adiffget(ad,mode)
% ADIFFGET returns parts of an adiff object. 
%
% x      = adiffget(A,'value')      returns the value of the adiff object A
% dx     = adiffget(A,'derivative') returns the jacobian
% [x,dx] = adiffget(A)              returns both the value and the jacobian
%
% The value of an adiff object is what value the object would have had it just
% been a vector. The jacobian is the derivative of the value with respect to the 
% variables implicitly created when the adiff object was first created.
%
% For example, a=adiff(1:10) creates an adiff object of 10 variables, with values 1 to 10.
% [x,dx] = adiffget(a) will then yield x=(1:10)', and dx is a 10 by 10 sparse identity matrix.
% If we then compute b = sum(a), b is an adiff object (because sum is defined for them).
% [x,dx] = adiffget(b) will yield x=55 i.e. sum(1:10), and dx =ones(1,10), since that is the 
% jacobian matrix for the sum.


if nargin==1
   x = ad.x;
   dx = ad.dx;
else
   switch mode(1)
   case 'v'
      x = ad.x;
   case 'd'
      x = ad.dx;
   end
end
