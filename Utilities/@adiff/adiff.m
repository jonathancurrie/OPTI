function ad = adiff(a,b,c)
% adiff(x) creates an autodiff object from a column vector x.
% adiff(x,'sparse') makes the adiff jacobian sparse
% adiff(x,dx,firstx) creates an adiff object from parts (internal use only)
%
% The adiff object that results has two parts, a value and a jacobian
% derivative. For example,
% a = adiff(1:10) yields a 10 dimensional adiff object. The value of
% a is 1:10, which are taken to be the values of 10 variables x(1)..x(10).
% The jacobian of a is the matrix dx(i)/dx(j), which is of course I.
%
% Following some calculations, we may have an adiff object b. The value of
% b is the result of the calculation, say some number 'y'. The jacobian of
% b is the row vector dy/dx(i), where x(i) are the original assumed variables. For
% example, if b=sum(a), then the value of b is 55, and the jacobian is ones(1,10).
%
% The value and jacobian can be retrieved using autoget.
%
% One very important restriction is that all the 
% calculations must be done starting with a single adiff object, otherwise
% the jacobian will not be computed correctly. For example, 
% a = adiff(1:10) yields a 10 dimensional adiff object
% b = adiff(1:3) yields a 3 dimensional adiff object; the three values
% of b are in no way related to the 10 values of a.


switch nargin
case 0
   z = struct('x',[],'dx',[],'root',[]);
case 1
   z = struct('x',a(:),'dx',eye(length(a)),'root',a(:));
case 2
   if ischar(b)&strcmp(b,'sparse')
      z = struct('x',a(:),'dx',speye(length(a)),'root',a(:));
   else
      error('Second argument must be ''sparse'' ');
   end
case 3
   z = struct('x',a(:),'dx',b,'root',c);
otherwise
   error('Wrong number of arguments');
end

ad = class(z,'adiff');