function [A,rl,ru] = gen2row(Ain,b,Aeq,beq)
%GEN2ROW Convert Linear Inequality & Equality Constraints to Linear A with Row Bounds
%   [A,rl,ru] = gen2row(Ain,b,Aeq,beq)

%   Copyright (C) 2011 Jonathan Currie (IPL)

if(isempty(Ain) && isempty(Aeq))
    A = []; rl = []; ru = [];
    return;
end
if(xor(isempty(Ain),isempty(b))), error('A + b cannot be empty!'); end
if(xor(isempty(Aeq),isempty(beq))), error('Aeq + beq cannot be empty!'); end

%Get sizes
nin = size(Ain,1);
neq = size(Aeq,1);
%Augment
A = [Ain;Aeq];
ru = [b;beq];
%Preallocate and index
rl = -Inf(nin+neq,1);
rl(nin+1:end) = beq;


