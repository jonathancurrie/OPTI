function [A,rl,ru,Aeq,beq] = rowe2gene(A,rl,ru,Aeq,beq)
%ROWE2GENE Move Linear Equalities in A with Row Bounds to Linear Equality Constraints
%   [A,rl,ru,Aeq,beq] = rowe2gene(A,rl,ru,Aeq,beq)

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(isempty(A))
    return;
end

%Get Indicies
eq = rl == ru;
if(~any(eq))
    return;
end
%Transpose as Required
if(size(rl,2) > 1)
    rl = rl';
end
if(size(ru,2) > 1)
    ru = ru';
end
if(size(beq,2) > 1)
    beq = beq';
end

%Create Equalities, augmenting to existing
Aeq = [Aeq; A(eq,:)]; beq = [beq; ru(eq)];

%Remove from row bounds
A =  A(~eq,:);
rl = rl(~eq);
ru = ru(~eq);
