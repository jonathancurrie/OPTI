function ff = fft(data)
% ff = fft(data) does the fourier transform of the data (power of 2 only)

% This m-file implementation of the fft is meant solely for use with adiff
% objects; it's far too inefficient for general use.

% recursive danielson-lanczos algorithm
n = length(data);
if n==1
   ff = data;
   return
end
ff0 = fft(subsref(data,struct('type','()','subs',{{1:2:n}})));
ff1 = fft(subsref(data,struct('type','()','subs',{{2:2:n}})));
w   = exp(-2*pi*i*(0:n/2-1)/n)';
ff  = [ff0+w.*ff1; ff0-w.*ff1];
