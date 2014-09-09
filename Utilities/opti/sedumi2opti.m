function prob = sedumi2opti(prob)
%SEDUMI2OPTI  Converts a SEDUMI Problem (At, b, C, K) to OPTI Format
% prob = sedumi2opti(prob)

%   Copyright (C) 2013 Jonathan Currie (I2C2)

if(~isfield(prob,'sdcone') || isempty(prob.sdcone) || ~isstruct(prob.sdcone))
    error('This function expects a prob.sdcone to contain a structure');
end
sdcone = prob.sdcone;
if((~isfield(sdcone,'A') && ~isfield(sdcone,'At')) || ~isfield(sdcone,'b') || (~isfield(sdcone,'C') && ~isfield(sdcone,'c')) || ~isfield(sdcone,'K'))
    error('This function expects prob.sdcone to contain the fields At (or A), b, C and K');
end
if(~isstruct(sdcone.K))
    error('SeDuMi prob.sdcone.K must be a structure!');
end
if(~isfield(sdcone.K,'s') && ~isfield(sdcone.K,'l'))
    error('SeDuMi prob.sdcone.K must contain the fields l or s (or both)');
end
if((isfield(sdcone.K,'scomplex') && ~isempty(sdcone.K.scomplex)) || (isfield(sdcone.K,'ycomplex') && ~isempty(sdcone.K.ycomplex)) || (isfield(sdcone.K,'xcomplex') && ~isempty(sdcone.K.xcomplex)))
    error('Complex variables are not supported in this interface - use SeDuMi directly to solve!'); 
end
if(isfield(sdcone.K,'f') && ~isempty(sdcone.K.f) && any(sdcone.K.f > 0))
    error('SeDuMi free variables are not supported in this interface - use SeDuMi directly to solve!'); 
end
if(isfield(sdcone.K,'q') && ~isempty(sdcone.K.q))
    error('SeDuMi SOCP constraints are not supported in this interface - use SeDuMi directly to solve!'); 
end
if(isfield(sdcone.K,'r') && ~isempty(sdcone.K.r))
    error('SeDuMi rotated cone constraints are not supported in this interface - use SeDuMi directly to solve!'); 
end

%Read SeDuMi problem, convert to OPTI standard primal form
% MINIMIZE   dot(f,x) SUCH THAT  SUM A_i*x_i - C >= 0. 
top = 1; dims = 0;
haveSDP = false; haveLP = false; 
K = sdcone.K;
if(isfield(sdcone,'C'))
    c = sdcone.C;
else
    c = sdcone.c;
end
if(isfield(sdcone,'A'))
    A = sdcone.A;
else
    A = sdcone.At;
end
b = full(sdcone.b);

%Check c, b is a vector
if(size(c,2) > 1 && size(c,1) > 1)
    error('SeDuMi c field must be a column vector!');
end
if(size(b,2) > 1 && size(b,1) > 1)
    error('SeDuMi b field must be a column vector!');
end

%Flip A based on b dimensions
if(size(A,2) ~= length(sdcone.b))
    if(size(A,1) == length(sdcone.b))
        A = A';
    else
        error('SeDuMi A (At) appears to  have incorrect dimensions.\nExpected one dimension of A to be the same as length as %s.','b');
    end
end
if(size(c,2) > size(c,1))
    c = c';
end
%Check length of C against A
if(length(c) ~= size(A,1))
    error('SeDuMi c does not have the same number of elements as rows in A (At)');
end
%Check dims entered in K
if(isfield(K,'s') && ~isempty(K.s) && K.s(1) > 0)
    haveSDP = true;
    dims = dims + sum(K.s.^2);
end
if(isfield(K,'l') && ~isempty(K.l) && K.l(1) > 0)
    if(length(K.l) > 1)
        error('The SeDuMi K.l field must contain a scalar, the number of linear constraints');
    end
    haveLP = true;
    dims = dims + K.l;
end
%Check dims against rows
if(dims ~= size(A,1))
    error('SeDuMi A (At) does not have correct dimensions as per the specifications in K. Expected %d rows.',dims);
end

%Extract Objective
prob.f = full(-b);

%Extract Linear Constraints
if(haveLP)
    n = K.l;
    if(isfield(prob,'A') && ~isempty(prob.A)) %concatenate
        prob.A = [prob.A; A(top:top+n-1,:)];
        prob.b = full([prob.b; prob.c(top:top+n-1,:)]);
    else
        prob.A = A(top:top+n-1,:);
        prob.b = full(c(top:top+n-1,:));
    end
    top = top+n;
end
%Extract Semidefinite Constraints
if (haveSDP) 
    prob.sdcone = cell(length(K.s),1);
    for i = 1:length(K.s)
        n = K.s(i);
        prob.sdcone{i} = [-c(top:top+n^2-1,:) -A(top:top+n^2-1,:)];
        top = top+n^2;
    end
else
    prob.sdcone = [];
end