function [Ain,bin,Aeq,beq] = mix2gen(A,b,e)
%MIX2GEN  Convert Mixed linear constraints to General Linear inequality and equality
%   [A,b,Aeq,beq] = mix2gen(A,b,e)
%
%   Use -1 for <=, 0 for =, and 1 for >= in vector e

%   Copyright (C) 2011 Jonathan Currie (I2C2)

if(size(A,1) ~= length(b))
    error('A and b sizes do not correspond');
end
if(length(b) ~= length(e))
    error('b and e are not the same length!');
end
if(any(e < -1) || any(e > 1))
    error('The vector e must only contain values -1, 0 or 1');
end

eq = find(e == 0);
leq = find(e == -1);
geq = find(e == 1);

%Transpose as neccesary
if(size(b,2) > 1)
    b = b';
end

%Process Equality Constraints
if(isempty(eq))
    Aeq = []; beq = [];
else
    Aeq = A(eq,:);
    beq = b(eq);
end

%Process <= Constraints
if(isempty(leq))
    Ain = []; bin = [];
else
    Ain = A(leq,:);
    bin = b(leq);
end

%Process >= Constraints
if(~isempty(geq))
    Ain = [Ain;-A(geq,:)];
    bin = [bin;-b(geq)];
end

