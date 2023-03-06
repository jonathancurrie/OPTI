function [ok,msg] = checkOptiSol(optObj,tol)
%CHECKOPTISOL  Check the solution to an optimization for errors
%
%   Called by OPTI / checkSol

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

ok = 1;
msg = sprintf('Solver Status:');

if(nargin < 2), tol = 1e-6; end

lintol = tol;
quadtol = tol;
nllintol = tol;
sdtol = tol;
inttol = 1e-5;

%Check Exit Flag
switch(optObj.ef)
    case 1
        msg = sprintf('%s Solved\n',msg);
    case 0
        msg = sprintf('%s Iterations / Time / Nodes Exceeded\n',msg);
        ok = 0;
    case -1
        msg = sprintf('%s Infeasible\n',msg);
        ok = 0;
    otherwise
        msg = sprintf('%s Solver Error\n',msg);
        ok = 0;
end

%Get Problem Properties
prob = optObj.prob; 
x = optObj.sol; [r,c] = size(x);
if(c > r), x = x'; end

%Convert SeDuMi problem to OPTI if specified
if(~isempty(prob.sdcone) && isstruct(prob.sdcone))
    prob = sedumi2opti(prob);
end

%Check Linear Constraints
if(~isempty(prob.rl))
    %Check Row Bounds
    if(~isempty(prob.A))
        cval = prob.A*x;       
        [msg,ok] = checkRowCon(cval,prob.rl,prob.ru,'Linear',lintol,msg,ok);        
    end
else
    %Check Inequality
    if(~isempty(prob.A))
        err = prob.A*x - prob.b;
        chk = err > lintol;
        if(any(chk))
            ok = 0;
            msg = sprintf('%s\nLinear <= Inequality Constraint(s) Broken: [tol %g]\n',msg,lintol);
            msg = solErrMsg(msg,true(size(chk)),chk,err);
        end
    end
    %Check Equality
    if(~isempty(prob.Aeq))
        err = abs(prob.Aeq*x - prob.beq);
        chk = err > lintol;
        if(any(chk))
            ok = 0;
            msg = sprintf('%s\nLinear == Equality Constraint(s) Broken: [tol %g]\n',msg,lintol);
            msg = solErrMsg(msg,true(size(chk)),chk,err);
        end
    end
end
%Check Bounds
if(~isempty(prob.lb))
    err = abs(x - prob.lb);
    chk = x < (prob.lb-lintol);
    if(any(chk))
        ok = 0;
        msg = sprintf('%s\nDecision Variable Lower Bound(s) Broken: [tol %g]\n',msg,lintol);
        msg = solErrMsg(msg,true(size(chk)),chk,err);
    end
end
if(~isempty(prob.ub))
    err = abs(x - prob.ub);
    chk = x > (prob.ub+lintol);
    if(any(chk))
        ok = 0;
        msg = sprintf('%s\nDecision Variable Upper Bound(s) Broken: [tol %g]\n',msg,lintol);
        msg = solErrMsg(msg,true(size(chk)),chk,err);
    end
end

%Check Quadratic Constraints
if(~isempty(prob.Q))
    if(iscell(prob.Q)) %multiple  
        %Evaluate each constraint
        cval = zeros(size(prob.qrl));
        for i = 1:length(prob.qrl)
           cval(i) = x'*prob.Q{i}*x + prob.l(:,i)'*x;            
        end
        %Check all constraints
        [msg,ok] = checkRowCon(cval,prob.qrl,prob.qru,'Quadratic',quadtol,msg,ok);
    else %single        
        cval = x'*prob.Q*x + prob.l'*x;      
        [msg,ok] = checkRowCon(cval,prob.qrl,prob.qru,'Quadratic',quadtol,msg,ok);
    end
end       

%Check Integer Constraints
if(any(prob.int.str == 'I'))
    idx = prob.int.str == 'I';
    err = abs(x - round(x));
    chk = err > inttol;
    if(any(chk(idx)))
        ok = 0;
        msg = sprintf('%s\nDecision Variable Integer Constraint(s) Broken: [tol %g]\n',msg,inttol);
        msg = solErrMsg(msg,idx,chk,err);
    end
end
%Check Binary Constraints
if(any(prob.int.str == 'B'))
    idx = prob.int.str == 'B';
    err = abs(max([x - round(x), x - 1, 0 - x],[],2));
    chk = err > inttol;
    if(any(chk(idx)))
        ok = 0;
        msg = sprintf('%s\nDecision Variable Binary Constraint(s) Broken: [tol %g]\n',msg,inttol);
        msg = solErrMsg(msg,idx,chk,err);
    end
end

%Check Semidefinite Constraints
if(~isempty(prob.sdcone))
    if(iscell(prob.sdcone))
        err = zeros(length(prob.sdcone),1);
        chk = zeros(length(prob.sdcone),1);
        for i = 1:length(prob.sdcone)
            cval = evalSDCone(prob.sdcone{i},x);
            err(i) = -min(eig(cval));
            chk(i) = err(i) > sdtol;            
        end
    else
        cval = evalSDCone(prob.sdcone,x);
        err = -min(eig(cval));
        chk = err > sdtol;
    end
    if(any(chk))
        ok = 0;
        msg = sprintf('%s\nSemidefinite Constraint(s) Broken: [tol %g]\n',msg,sdtol);
        msg = solErrMsg(msg,true(size(chk)),chk,err);
    end
end

%Check Nonlinear Constraints
if(~isempty(prob.nlcon))
    %Convert general constraints to row (if not already)
    if(~isempty(prob.nlrhs)), prob = nmix2row(prob); end
    %Get cval
    cval = prob.nlcon(x);
    %Run check
    [msg,ok] = checkRowCon(cval,prob.cl,prob.cu,'Nonlinear',nllintol,msg,ok);
end

function [msg,ok] = checkRowCon(cval,rl,ru,str,tol,msg,ok)
%Ensure columns
cval = cval(:);
rl = rl(:);
ru = ru(:);
%Indices
eq = rl == ru; neq = ~eq;
ile = isfinite(ru) & neq;
ige = isfinite(rl) & neq;
if(any(ile))
    err = (cval - ru);
    chk = err > tol;
    if(any(chk(ile)))
        ok = 0;
        msg = sprintf('%s\n%s <= Inequality Constraint(s) Broken: [tol %g]\n',msg,str,tol);
        msg = solErrMsg(msg,ile,chk,err);
    end
end
if(any(ige))
    err = (rl - cval);
    chk = err > tol;
    if(any(chk(ige)))
        ok = 0;
        msg = sprintf('%s\n%s >= Inequality Constraint(s) Broken: [tol %g]\n',msg,str,tol);
        msg = solErrMsg(msg,ige,chk,err);
    end
end
if(any(eq))
    err = abs(cval - rl);
    chk = err > tol;
    if(any(chk(eq)))
        ok = 0;
        msg = sprintf('%s\n%s == Equality Constraint(s) Broken: [tol %g]\n',msg,str,tol);
        msg = solErrMsg(msg,eq,chk,err);
    end
end

function msg = solErrMsg(msg,indx,chk,err)
for i = 1:length(chk)
    if(indx(i) == true && chk(i) == true)
        msg = sprintf('%s #%-5d Error: %-10g\n',msg,i,err(i));
    end
end

function cval = evalSDCone(cone,x)
dim = sqrt(size(cone,1));
cval = -cone(:,1);
for i = 1:length(x)
    cval = cval + x(i).*cone(:,i+1);
end
cval = reshape(cval,dim,dim);
