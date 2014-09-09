function [x,fval,flag,info] = mosekRes(prob,rcode,res,t)
% MOSEKRES  Extract solution and status from MOSEK results
%
%   [x,fval,flag,info] = mosekRes(prob,rcode,res,t)

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (I2C2)

hasSol = isfield(res,'sol');
if(hasSol && isfield(res.sol,'itr'))
    isint = 0;
else
    isint = 1;
end

%Get Solution
if(hasSol)
    if(isint)
        x = res.sol.int.xx;
        fval = res.sol.int.pobjval;
    else
        x = res.sol.itr.xx;
        fval = res.sol.itr.pobjval;
    end
else
    x = [];
    fval = [];
end

%Get symbolic constants
if(isfield(res,'symbcon'))
    sc = res.symbcon;
else
    [~,r] = mosekopt('symbcon');
    sc = r.symbcon;
end

%Check for convexity or license errors
switch(rcode)
    case sc.MSK_RES_ERR_INV_PROBLEM
        flag = -3;
        info = struct('Iterations',[],'Time',toc(t),'Algorithm','MOSEK','Status','Invalid Problem - Possibly nonconvex?');
        return;
    case {sc.MSK_RES_ERR_MISSING_LICENSE_FILE, sc.MSK_RES_ERR_LICENSE,sc.MSK_RES_ERR_LICENSE_EXPIRED,sc.MSK_RES_ERR_LICENSE_VERSION, ...
          sc.MSK_RES_ERR_SIZE_LICENSE ,sc.MSK_RES_ERR_PROB_LICENSE,sc.MSK_RES_ERR_FILE_LICENSE,sc.MSK_RES_ERR_MISSING_LICENSE_FILE, ...
          sc.MSK_RES_ERR_SIZE_LICENSE_CON,sc.MSK_RES_ERR_SIZE_LICENSE_VAR,sc.MSK_RES_ERR_SIZE_LICENSE_INTVAR,sc.MSK_RES_ERR_OPTIMIZER_LICENSE, ...
          sc.MSK_RES_ERR_FLEXLM,sc.MSK_RES_ERR_LICENSE_SERVER,sc.MSK_RES_ERR_LICENSE_MAX,sc.MSK_RES_ERR_LICENSE_FEATURE,sc.MSK_RES_ERR_LICENSE_CANNOT_CONNECT, ...
          sc.MSK_RES_ERR_LICENSE_INVALID_HOSTID,sc.MSK_RES_ERR_LICENSE_SERVER_VERSION,sc.MSK_RES_ERR_OPEN_DL ,sc.MSK_RES_ERR_OLDER_DLL,sc.MSK_RES_ERR_NEWER_DLL}
        flag = -4;
        info = struct('Iterations',[],'Time',toc(t),'Algorithm','MOSEK','Status',sprintf('License / Program Error: %s',res.rmsg));
        return;
end

%Otherwise extract exit information
if(hasSol)
    if(isint)
        solsta = res.sol.int.solsta; %solution status
        prosta = res.sol.int.prosta; %problem status
    else
        solsta = res.sol.itr.solsta;
        prosta = res.sol.itr.prosta;
    end
end %no else required as below checks as well   

dual = 0; %unsure why this is needed

if(hasSol)
    if(rcode == 0 && solsta == sc.MSK_SOL_STA_OPTIMAL)
        flag = 1;
        status = 'Optimal';
    elseif(rcode == 0 && solsta == sc.MSK_SOL_STA_INTEGER_OPTIMAL)
        flag = 1;
        status = 'Integer Optimal';
    elseif(~dual && prosta == sc.MSK_PRO_STA_PRIM_INFEAS || dual && prosta == sc.MSK_PRO_STA_PRIM_AND_DUAL_INFEAS)
        flag = -1;
        status = 'Primal Infeasible';
    elseif(~dual && prosta == sc.MSK_PRO_STA_DUAL_INFEAS || dual && prosta == sc.MSK_PRO_STA_PRIM_INFEAS)
        flag = -2;
        status = 'Dual Infeasible';    
    else
        flag = -5;
        if(~isempty(res.rmsg))
            status = sprintf('MOSEK Status: %s',res.rmsg);
        elseif(isfield(res,'rcodestr'))
            status = sprintf('MOSEK Status: %s',res.rcodestr);
        else
            status = 'MOSEK Status: Unknown (Enable display to examine problem)';
        end
    end    
elseif(rcode == sc.MSK_RES_TRM_MAX_ITERATIONS)
    flag = 0;
    status = 'Exceeded Iterations'; 
else
    flag = -6;
    status = sprintf('MOSEK Error: %s',res.rmsg); %not very informative...
end    

%Build info structure
if(isfield(res,'info'))   
    if(isfield(res,'sol') && isfield(res.sol,'itr'))
        info.Iterations = res.info.MSK_IINF_INTPNT_ITER;
    else
        info.Iterations = res.info.MSK_IINF_MIO_NUM_BRANCH;
    end
else
    info.Iterations = [];
end
info.Time = toc(t);
info.Algorithm = 'MOSEK';

%Assign Status
info.Status = status;

%Assign lagrangian at solution
nineq = size(prob.a,1); %includes quadratics
if(isfield(res,'sol') && isfield(res.sol,'itr')) %MI problems do not have dual information        
    info.Lambda.lower = res.sol.itr.slx;
    info.Lambda.upper = res.sol.itr.sux;
    if(isfield(res.sol.itr,'barx') && isfield(prob,'bardim'))
        info.DualObjective = -fval;
        fval = -fval; %return to primal form        
        x = -res.sol.itr.y; %saved in dual solution
        ndim = length(prob.bardim);
        info.Lambda.X = cell(ndim,1); 
        info.Lambda.S = cell(ndim,1); 
        top = 1; barx = res.sol.itr.barx; bars = res.sol.itr.bars;
        for i = 1:ndim
            n = prob.bardim(i);
            info.Lambda.X{i} = dmat(barx(top:top+n*(n+1)/2-1));
            info.Lambda.S{i} = dmat(bars(top:top+n*(n+1)/2-1));
            top = top + n*(n+1)/2;
        end
    else
        info.Lambda.ineqlin = -res.sol.itr.y(1:nineq);
        info.Lambda.eqlin = -res.sol.itr.y(nineq+1:end);
    end
else
    info.Lambda = [];
end
    


function [A]=dmat(V)
% DSDP5.0 
% Copyright (c) 2003 by
% S. Benson and Y. Ye
% Last modified: December 2003
[m,m1]=size(V);
if (m1~=1), error('Not a column vector.'); end;
n=floor(sqrt(2*m));
if (n*(n+1)/2 ~= m) 
     error('Impossible Dimension. Cannot convert to square matrix.'); 
end;
nzV=nnz(V);
if (issparse(V) || nzV < n*(n+1)/2), A=sparse(n,n); else A=zeros(n,n);  end;
for i=1:n, 
    k1=i*(i-1)/2+1;  k2=i*(i+1)/2;
    A(1:i,i)= V(k1:k2); 
end;
A=A+triu(A,1)';    
    