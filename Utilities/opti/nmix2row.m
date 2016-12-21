function prob = nmix2row(prob)
%NMIX2ROW  Convert Mixed Nonlinear Constraints to Nonlinear Row bounds
%   prob = nmix2row(prob)
%
%   Use -1 for <=, 0 for =, and 1 for >= in vector e

%   Copyright (C) 2012 Jonathan Currie (IPL)

%Assign common vars
rhs = prob.nlrhs;
e = prob.nle;

%Transpose as Required
if(size(rhs,2) > 1)
    rhs = rhs';
end
if(size(e,2) > 1)
    e = e';
end
if(xor(isempty(rhs),isempty(e)))
    error('You must supply both nlrhs and nle!');
end
if(length(rhs) ~= length(e))
    error('nlrhs and nle are not the same length!');
end
if(any(e < -1) || any(e > 1))
    error('The vector e must only contain values -1, 0 or 1');
end

%Defaults
ncon = length(e);
prob.cl = -Inf(ncon,1);
prob.cu = Inf(ncon,1);

%Indices
ieq = e == 0;
icl = (e ==  1 | ieq);
icu = (e == -1 | ieq);

%Fill In
prob.cl(icl) = rhs(icl);
prob.cu(icu) = rhs(icu);

%Remove unused fields
prob.nlrhs = []; prob.nle = [];  

