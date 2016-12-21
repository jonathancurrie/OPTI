function [A,b,Aeq,beq] = row2gen(Ain,rl,ru)
%ROW2GEN Convert Linear A with Row Bounds to Linear Inequality & Equality Constraints
%   [A,b,Aeq,beq] = row2gen(A,rl,ru)

%   Copyright (C) 2012 Jonathan Currie (IPL)

if(isempty(Ain))
    A = []; b = []; Aeq = []; beq = [];
    return;
elseif(isempty(ru))
    error('You must supply Ain, rl and ru!');
end
if(isempty(rl))
    rl = -Inf(size(ru)); %default
end
%Transpose as Required
if(size(rl,2) > 1)
    rl = rl';
end
if(size(ru,2) > 1)
    ru = ru';
end

%Indices
eq = rl == ru; neq = ~eq;
ile = isfinite(ru) & neq;
ige = isfinite(rl) & neq;

%Ineq
A = [ Ain(ile,:);
     -Ain(ige,:)];
b = [ ru(ile);
     -rl(ige)];
%Eq
if(any(eq))
    Aeq = Ain(eq,:); beq = ru(eq);
else
    Aeq = []; beq = [];
end

