function [v,vineq,veq] = optiConViolation(x,A,b,Aeq,beq,nlcon,nlrhs,nle,Q,l,qrl,qru,lb,ub,int,sdcone)
%OPTICONVIOLATION  Evaluate ALL OPTI constraints and return absolute violation

v = []; vineq = []; veq = [];

%Bounds
if(~isempty(lb))
    lb = lb(~isinf(lb));
    vl = x - lb; vl(vl > 0) = 0;
    v = [v;vl];
end
if(~isempty(ub))
    ub = ub(~isinf(ub));
    vu = x - ub; vu(vu < 0) = 0;
    v = [v;vu];
end
%Linear Constraints
if(~isempty(A))
    vl = A*x - b; vl(vl < 0) = 0;
    v = [v;vl];
    vineq = [vineq;vl];
end
if(~isempty(Aeq))
    v = [v;Aeq*x - beq];
    veq = [veq;v];
end
%Quadratic Constraints
if(~isempty(Q))
    q = x'*Q*x + l'*x;
    idxe = qrl == qru;    
    idxl = ~isinf(qrl) & ~idxe;
    idxu = ~isinf(qru) & ~idxe;
    vl = q(idxl) - qrl(idxl); vl(vl > 0) = 0;
    vu = q(idxu) - qru(idxu); vu(vu < 0) = 0;
    ve = q(idxe) - qrl(idxe);
    v = [v;vl;vu];
    vineq = [vineq;vl;vu];
    veq = [veq;ve];
end
%Nonlinear Constraints
if(~isempty(nlcon))
    n = nlcon(x);
    ile = nle == -1;
    ige = nle == 1;
    ieq = nle == 0;
    cel = n(ile); vl = cel(:) - nlrhs(ile); vl(vl < 0) = 0;
    ceu = n(ige); vu = ceu(:) - nlrhs(ige); vu(vu > 0) = 0;
    cee = n(ieq); ve = cee(:) - nlrhs(ieq);
    v = [v;vl;vu;ve];
    vineq = [vineq;vl;vu];
    veq = [veq;ve];
end
%Semidefinite Constraints
if(~isempty(sdcone))
    error('not implemented');
end

%Integer Constraints


%Absolute Value for Violation
v = abs(v);
vineq = abs(vineq);
veq = abs(veq);

