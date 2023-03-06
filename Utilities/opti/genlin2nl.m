function prob = genlin2nl(prob,sparsity,warn)
%GENLIN2NL Convert General Linear Constraints to Mixed Nonlinear Constraints

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

if(isempty(prob.A) && isempty(prob.Aeq))
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
    if(isempty(prob.Aeq))        
        prob.nlcon = @(x) prob.A*x;
        if(sparsity)
            prob.nljac = @(x) sparse(prob.A);
            prob.nljacstr = @() sparse(double(prob.A ~= 0));
        else
            prob.nljac = @(x) prob.A;
            prob.nljacstr = @() double(prob.A ~= 0);
        end
    elseif(isempty(prob.A))
        prob.nlcon = @(x) prob.Aeq*x;
        if(sparsity)
            prob.nljac = @(x) sparse(prob.Aeq);
            prob.nljacstr = @() sparse(double(prob.Aeq ~= 0));
        else
            prob.nljac = @(x) prob.Aeq;
            prob.nljacstr = @() double(prob.Aeq ~= 0);
        end
    else
        prob.nlcon = @(x) [prob.A*x;
                           prob.Aeq*x];
        if(sparsity)
            prob.nljac = @(x) sparse([prob.A; prob.Aeq]);
            prob.nljacstr = @() sparse([double(prob.A ~= 0); double(prob.Aeq ~= 0)]);  
        else
            prob.nljac = @(x) [prob.A; prob.Aeq];
            prob.nljacstr = @() [double(prob.A ~= 0); double(prob.Aeq ~= 0)];  
        end
    end
    prob.nlrhs = [prob.b;prob.beq];
    prob.nle = [-1*ones(size(prob.b));zeros(size(prob.beq))];
    %Add sizes       
    prob.sizes.nnlineq = length(prob.b);
    prob.sizes.nnleq = length(prob.beq); 
    
    prob.lincon = 1; %indicate all linear constraints

%Otherwise Augment linear to nonlinear constraints
else   
    %Augment nonlinear constraint function + jacobian + jac structure
    if(isempty(prob.Aeq))
        prob.nlcon = @(x) [prob.nlcon(x);
                           prob.A*x];
        prob.nljac = @(x) [prob.nljac(x);
                           prob.A];
        if(~isempty(prob.nljacstr))
            prob.nljacstr = @() [prob.nljacstr(); sparse(double(prob.A ~= 0))];
        end
    elseif(isempty(prob.A))
        prob.nlcon = @(x) [prob.nlcon(x);
                           prob.Aeq*x];
        prob.nljac = @(x) [prob.nljac(x);
                           prob.Aeq];
        if(~isempty(prob.nljacstr))
            prob.nljacstr = @() [prob.nljacstr(); sparse(double(prob.Aeq ~= 0))];
        end
    else
        prob.nlcon = @(x) [prob.nlcon(x);
                           prob.A*x;
                           prob.Aeq*x];
        prob.nljac = @(x) [prob.nljac(x);
                           prob.A;
                           prob.Aeq];
        if(~isempty(prob.nljacstr))
            prob.nljacstr = @() [prob.nljacstr(); sparse(double(prob.A ~= 0)); sparse(double(prob.Aeq ~= 0))];
        end
    end
    %Augment rhs + type
    prob.nlrhs = [prob.nlrhs;
                  prob.b;
                  prob.beq];
    prob.nle = [prob.nle;
                -1*ones(size(prob.b));
                zeros(size(prob.beq))];
    %Add sizes       
    prob.sizes.nnlineq = prob.sizes.nnlineq + length(prob.b);
    prob.sizes.nnleq = prob.sizes.nnleq + length(prob.beq);                    

    prob.lincon = 0; %not all linear    
end

if(warn > 1)
    optiwarn('opti:lin_nl','General Linear Constraints have been converted to Mixed Nonlinear Constraints to suit NLP Solver interface');
end

end

