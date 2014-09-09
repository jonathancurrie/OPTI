function ff = ifft(data)
% ff = ifft(data) does the inverse fourier transform of the data (power of 2 only)

% This m-file implementation of the fft is meant solely for use with adiff
% objects; it's far too inefficient for general use.

ff = fft(conj(data))/length(data);
