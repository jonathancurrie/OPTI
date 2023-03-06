function [At,b,c,K] = opti2sedumi(prob)
%OPTI2SEDUMI  Converts an OPTI Problem (f,A,b,sdcone) to SeDuMi Format
% [At,b,c,K] = opti2sedumi(prob)

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

%Check not already in sedumi format
if(isfield(prob,'sdcone') && ~isempty(prob.sdcone) && isstruct(prob.sdcone))
    At = getSDMIField(prob.sdcone,'At','A');
    b  = getSDMIField(prob.sdcone,'b','B');
    c  = getSDMIField(prob.sdcone,'c','C');
    K  = getSDMIField(prob.sdcone,'K','k');
    return;
end

%Check we have a linear f
if(~isfield(prob,'f') || isempty(prob.f) || ~isnumeric(prob.f))
    error('The SDP to convert must have a linear objective f in prob.f');
end
b = -prob.f; ndec = length(b);

%Convert linear constraints to general format
if(~isempty(prob.rl) || ~isempty(prob.ru))
    [AA,bb,Aeq,beq] = row2gen(prob.A,prob.rl,prob.ru);
else
    AA = prob.A; bb = prob.b; Aeq = prob.Aeq; beq = prob.beq;
end
%Convert all constraints and bounds to linear inequalities
[At,c] = inequal(AA,bb,Aeq,beq,prob.lb,prob.ub);
K.l = size(At,1);

%Check for Semidefinite constraints, concat as required
if(~isempty(prob.sdcone))
    if(iscell(prob.sdcone))
        for i = 1:length(prob.sdcone)
            chkSDDim(prob.sdcone{i},ndec,i);
            K.s(i) = sqrt(size(prob.sdcone{i},1));
            c = [c;-prob.sdcone{i}(:,1)]; %#ok<AGROW>
            At = [At;-prob.sdcone{i}(:,2:end)]; %#ok<AGROW>
        end
    else
        chkSDDim(prob.sdcone,ndec,1);
        K.s = sqrt(size(prob.sdcone,1));
        c = [c;-prob.sdcone(:,1)];
        At = [At;-prob.sdcone(:,2:end)];
    end
end

%Check for quadratic constraints
if(~isempty(prob.Q)), optiwarn('opti:sedumiQC','Currently quadratic constraints are not converted to SeDuMi format'); end


function chkSDDim(cone,ndec,i)
if(size(cone,2) ~= ndec+1)
    error('Semidefinite Cone %d does not have the correct number of columns',i);
end
dim = sqrt(size(cone,1));
if(floor(dim) ~= dim)
    error('Semidefinite Cone %d does not have the correct number of rows to form a square matrix');
end


function val = getSDMIField(conestr,f1,f2)
if(isfield(conestr,f1))
    val = conestr.(f1);
elseif(isfield(prob.sdcone,f2))
    val = conestr.(f2);
else
    error('It appears the problem is already in SeDuMi format, but OPTI could not find prob.sdcone.%s',f1);
end