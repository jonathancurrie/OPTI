function c = mpower(a,b)
% MPOWER implements c = a^b, where either a or b is an adiff object.
% For adiff objects, a^b is taken to be equivalent to a.^b

c=power(a,b);
