function prob = rowlin2nl(prob,sparsity,warn)
%ROWLIN2NL Convert Row Based Linear Constraints to Row Based Nonlinear Constraints

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(isempty(prob.A))
    return; %nothing to do
end

if(~exist('warn','var') || isempty(warn))
    warn = 1;
end
if(~exist('sparsity','var') || isempty(sparsity))
    sparsity = 1;
end

%If we don't have any nonlinear constraints, just build from linear ones
if(isempty(prob.nlcon))
    %Build nonlinear constraints + jacobian + jac structure (assumes jac + str not provided)        
    prob.nlcon = @(x) prob.A*x;
    if(sparsity)
        prob.nljac = @(x) sparse(prob.A);
        prob.nljacstr = @() sparse(double(prob.A ~= 0));
    else
        prob.nljac = @(x) prob.A;
        prob.nljacstr = @() double(prob.A ~= 0);
    end
    %Assign constraint bounds
    prob.cl = prob.rl;
    prob.cu = prob.ru;
    %Add sizes       
    eq = prob.cl == prob.cu; neq = ~eq;
    ile = isfinite(prob.cu) & neq;
    ige = isfinite(prob.cl) & neq;
    prob.sizes.nnlineq = sum(ile + ige);
    prob.sizes.nnleq = sum(eq);
    %Indicate all linear
    prob.lincon = 1;

%Otherwise Augment linear to nonlinear constraints
else   
    %Augment nonlinear constraint function   
    prob.nlcon = @(x) [prob.nlcon(x);
                       prob.A*x];
    %Augment nonlinear Jacobian if it exists
    if(~isempty(prob.nljac))
        if(sparsity)
            A = sparse(prob.A);
            prob.nljac = @(x) [prob.nljac(x);
                               A];
        else
            prob.nljac = @(x) [prob.nljac(x);
                               prob.A];
        end
    end
    %Augment nonlinear jacobian structure if it exists
    if(~isempty(prob.nljacstr))
        if(sparsity)
            A = sparse(prob.A ~= 0);
            prob.nljacstr = @() [prob.nljacstr(); A];
        else
            prob.nljacstr = @() [prob.nljacstr(); double(prob.A ~= 0)];
        end
    end
    
    %Augment rhs + type
    prob.cl = [prob.cl;
               prob.rl];
    prob.cu = [prob.cu;
               prob.ru];
    %Add sizes       
    eq = prob.rl == prob.ru; neq = ~eq;
    ile = isfinite(prob.ru) & neq;
    ige = isfinite(prob.rl) & neq;
    prob.sizes.nnlineq = prob.sizes.nnlineq + sum(ile + ige);
    prob.sizes.nnleq = prob.sizes.nnleq + sum(eq);                  
    %Indicate not all linear
    prob.lincon = 0;    
end

if(warn > 1)
    optiwarn('opti:rowlin_nl','RowLin2NL - Linear Constraints have been converted to Nonlinear Constraints to suit NLP Solver interface');
end

end

