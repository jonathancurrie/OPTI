function [A,b] = inequal(Ain,bin,Aeq,beq,lb,ub)
%INEQUAL Convert Equality and Bounds to Inequality Constraints
%   [A,b] = inequal(A,b,Aeq,beq,lb,ub)

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

if(isempty(Ain)), Ain = []; end
if(isempty(bin)), bin = []; end

%Process Equality Constraints
if(isempty(Aeq))
    Ae = []; be = [];
else
    Ae = [Aeq;-Aeq];
    be = [beq;-beq];
end

if(nargin > 4)
    %Process Bounds
    Ab = -1*eye(length(lb));
    Ab = [Ab; eye(length(ub))];

    if(size(lb,2) > 1)
        lb = lb';
    end
    if(size(ub,2) > 1)
        ub = ub';
    end
    bb = [-lb;ub];

    %Remove infinite bounds
    ind = find(isinf(bb));
    Ab(ind,:) = [];
    bb(ind) = [];
    %Augment Bounds + Equality
    A = [Ain;Ae;Ab]; b = [bin;be;bb];
else
    A = [Ain;Ae]; b = [bin;be];
end



