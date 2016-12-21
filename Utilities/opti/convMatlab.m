function mprob = convMatlab(prob,opts)
%CONVMATLAB Convert OPTI problem to MATLAB Optimization Toolbox problem
%
%   mprob = convMatlab(prob,opts)

%   Copyright (C) 2011 Jonathan Currie (IPL)

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Make sure we have a suitable problem
if(~any(strcmpi(prob.type,{'NLP','UNO','SNLE','NLS','DNLS'})))
    error('You can only convert NLSs, UNOs, and NLPs to MATLAB format!');
end

%Common Args
mprob.options = opts.solverOpts;

%Check for iterfun
if(isfield(opts,'iterfun') && ~isempty(opts.iterfun))
    mprob.options.OutputFcn = @(x,o,s) mCall(x,o,s,opts.iterfun);
end

%Starting Guess
if(isfield(prob,'x0'))
    mprob.x0 = prob.x0;
end

%Problem dependent args
switch(upper(prob.type))
    case 'UNO'
        mprob.solver = 'fminsearch';
        mprob.objective = prob.fun;
    
    case 'SNLE'
        mprob.solver = 'fsolve';
        mprob.options = optimset(mprob.options);
        %Check for Gradient
        if(~isempty(prob.f))
            %Set gradient function
            mprob.objective = @(x) gradObj(x,prob.fun,prob.f);                
            mprob.options = optimset(mprob.options,'Jacobian','on');                                  
        else %skip grad
            mprob.objective = prob.fun; 
        end
        
    case {'DNLS','NLS'}
        mprob.solver = 'lsqnonlin';
        mprob.options = optimset(mprob.options);
        
        %Bounds stuff
        mprob.lb = prob.lb;
        mprob.ub = prob.ub;
        %Gen Objective Function
        if(~isempty(prob.xdata) && nargin(prob.fun)==2) %assume xdata supplied as 2nd arg to fun
            if(~isempty(prob.ydata))
                fun = @(x) prob.fun(x,prob.xdata)-prob.ydata;
            else
                fun = @(x) prob.fun(x,prob.xdata);
            end
        else
            if(~isempty(prob.ydata))
                fun = @(x) prob.fun(x)-prob.ydata;
            else
                fun = @(x) prob.fun(x);
            end
        end
        prob.fun = fun;
        %Check for Gradient
        if(~isempty(prob.f))
            if(~isempty(strfind(char(prob.f),'mklJac')) && ~isempty(strfind(char(prob.f),'prob.fun')))
                %Assume missing xdata input, fix here
                grad = @(x)mklJac(fun,x,length(prob.ydata));
            else
                if(nargin(prob.f) == 2) %assume also requires xdata
                    grad = @(x) prob.f(x,prob.xdata);
                else
                    grad = prob.f;
                end
            end
            %Set gradient function
            mprob.objective = @(x) gradObj(x,fun,grad);                
            mprob.options = optimset(mprob.options,'Jacobian','on');                                  
        else %skip grad
            mprob.objective = prob.fun; 
        end
        
    case 'NLP'       
        mprob.solver = 'fmincon';       
        %IP is the only one that allows all options
        if(isempty(mprob.options) || ~isfield(mprob.options,'Algorithm') || isempty(mprob.options.Algorithm))
            mprob.options = optimset(mprob.options,'Algorithm','interior-point'); 
        end
        mprob.algorithm = mprob.options.Algorithm;                        
        %Linear stuff
        mprob.Aineq = prob.A;
        mprob.bineq = prob.b;
        mprob.Aeq = prob.Aeq;
        mprob.beq = prob.beq;
        mprob.lb = prob.lb;
        mprob.ub = prob.ub;
        %Check for Gradient
        if(~isempty(prob.f))
            %Set gradient function
            mprob.objective = @(x) gradObj(x,prob.fun,prob.f);                
            mprob.options = optimset(mprob.options,'GradObj','on');                                               
        else %skip grad
            mprob.objective = prob.fun; 
        end
        
        %Check for NL Constraints
        if(~isempty(prob.nlcon))
            %Get Constraint Types
            max_in = prob.nle == 1;
            min_in = prob.nle == -1;
            eq = prob.nle == 0;
            nc = length(prob.nle);
            %Check for NL Jac            
            if(~isempty(prob.nljac)) 
                mprob.nonlcon = @(x) gradNlCon(x,prob.nlcon,prob.nljac,prob.nlrhs,max_in,min_in,eq);
                mprob.options = optimset(mprob.options,'GradConstr','on');  
            else %skip jac              
                mprob.nonlcon = @(x) nlCon(x,prob.nlcon,prob.nlrhs,max_in,min_in,eq);
            end
        else
            max_in = []; min_in = []; eq = []; nc = 0; %used below
        end
        
        %Check for Hessian
        if(~isempty(prob.H))
            %Determine type of Hessian passed (sym/tril/triu)
            hfmt = optiMatType(optiGenTestHessian(prob));
            H = @(x,lambda) mHess(prob.H,x,lambda,max_in,min_in,eq,nc,hfmt);
            %Set Hessian        
            mprob.options = optimset(mprob.options,'Hessian','user-supplied','HessFcn',H);
        elseif(strcmpi(mprob.algorithm,'trust-region-reflective'))
            mprob.options = optimset(mprob.options,'Hessian','fin-diff-grads'); 
        else
            mprob.options = optimset(mprob.options,'Hessian','lbfgs'); %limited memory BFGS update (10 iter)
        end
        
    otherwise
        error('Not implemented yet');
end



function [j,dj] = gradObj(x,fun,grad)
% Handle to allow matlab to get objective and gradient in a single function
j = fun(x);
if(nargout > 1)
    dj = grad(x);    
end    

function [cin,ceq] = nlCon(x,fun,rhs,max_in,min_in,eq)
% Handle to allow matlab to get nonlinear inequality and equality
% constraints in a single function with selectable bounds 
% (not very efficient)

%Get Constraint Eval
sol = fun(x);
%Defaults
cin = [];
ceq = [];
%Assign results with bounds
if(any(max_in))
    cin = -sol(max_in) + rhs(max_in);
end
if(any(min_in))
    cin = [cin; sol(min_in) - rhs(min_in)];
end
if(any(eq))
    ceq = sol(eq) - rhs(eq);
end

function [cin,ceq,dcin,dceq] = gradNlCon(x,fun,jac,rhs,max_in,min_in,eq)
% Handle to allow constraint and jacobian for nonlinear constraints

%Get Fvals
[cin,ceq] = nlCon(x,fun,rhs,max_in,min_in,eq);
%If gradient required, calculate
if(nargout > 2)
    sol = jac(x);
    %Defaults
    dcin = [];
    dceq = [];
    %Assign results with correct sign
    if(any(max_in))
        dcin = -sol(max_in,:);
    end
    if(any(min_in))
        dcin = [dcin; sol(min_in,:)];
    end
    if(any(eq))
        dceq = sol(eq,:);
    end
    %MATLAB expects transposed Jac
    dcin = dcin.';
    dceq = dceq';
end

function H = mHess(fun,x,lambda,max_in,min_in,eq,nc,hfmt)
% Handle to convert Matlab Hessian to OPTI Hessian

%Get Lambda as an array (assume Matlab has not reordered constraints)
l = zeros(nc,1);
if(nc)
   l(eq) = lambda.eqnonlin; 
   l((max_in | min_in)) = lambda.ineqnonlin; 
   l(max_in) = -l(max_in); %flip for >= constraints
end

%Evaluate callback
switch(nargin(fun))
    case 1
        H = fun(x);
    case 2
        H = fun(x,l);
    case 3
        H = fun(x,1,l);
    otherwise
        error('Unknown Hessian callback format, only 1-3 arguments are supported');
end
%Adjust as Matlab requires a full symmetric Hessian
switch(lower(hfmt))
    case 'l'
        H = H + tril(H,-1).';
    case 'u'
        H = H + triu(H,1).';
end


function stop = mCall(x,optimValues,state,fun)
% Handle to convert MATLAB callback to OPTI callback

stop = false;
switch state
    case 'iter'
        if isfield(optimValues,'fval')
            stop = fun(optimValues.iteration,optimValues.fval,x);
        else
            if(isscalar(optimValues.residual))
                stop = fun(optimValues.iteration,optimValues.residual,x);
            else
                stop = fun(optimValues.iteration,optimValues.resnorm,x);
            end
        end
        drawnow;
end



