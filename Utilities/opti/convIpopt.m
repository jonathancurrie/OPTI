function nlprob = convIpopt(prob,opts)
%CONVIPOPT Convert OPTI problem to IPOPT Nonlinear problem
%
%   nlprob = convIpopt(prob,opts)

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Get Warning Level
if(strcmpi(opts.warnings,'all'))
    warn = 2;
elseif(strcmpi(opts.warnings,'critical'))
    warn = 1;
else
    warn = 0;
end

%Check if the user wants to use IPOPT QP/QCQP solver
if(~isempty(strfind(prob.type,'QP')))
    if(warn > 1)
        optiwarn('OPTI:QCQP2NLP','Quadratic Program converted to (MI)NLP to suit IPOPT/BONMIN interface');
    end
    if(prob.sizes.nqc==0) %remember hessian only constant for LP/QP + mehrotra can crash on QCQP
        opts.solverOpts = ipoptset(opts.solverOpts,'hessian_constant','yes','mehrotra_algorithm','yes','mu_oracle','probing');
    end    
    [prob,opts] = QCQP2NLP(prob,opts);    
end

%Fillers
if(~isempty(prob.rl)), linvec = false(length(prob.rl),1); else linvec = logical([]); end
if(~isempty(prob.cl)), nlvec = false(length(prob.cl),1); else nlvec = logical([]); end

%Generate Nonlinear Constraint Indices
if(~isempty(prob.cl))
    eq = prob.cl == prob.cu; neq = ~eq;
    ile = isfinite(prob.cu) & neq;
    ige = isfinite(prob.cl) & neq;    
    options.nleq = [eq; linvec];
    options.nlineq = [ile | ige; linvec];
else
    options.nleq = [];
    options.nlineq = [];
end
%Linear Constraints will be Augmented, check and add indices
if(~isempty(prob.rl))
    eq = prob.rl == prob.ru; neq = ~eq;
    ile = isfinite(prob.ru) & neq;
    ige = isfinite(prob.rl) & neq;    
    options.eq = [nlvec; eq];
    options.ineq = [nlvec; ile | ige];
else
    options.eq = [];
    options.ineq = [];
end
 
%Configure Problem
if(~isempty(prob.lb)), options.lb = prob.lb; end
if(~isempty(prob.ub)), options.ub = prob.ub; end
if(~isempty(prob.cl)), options.cl = prob.cl; end
if(~isempty(prob.cu)), options.cu = prob.cu; end

% Set IPOPT options (optiset options override common ones)
%Build ipoptset struct, adding any passed options (ignored by ipoptset if not recognised)
options.ipopt = ipoptset(opts.solverOpts);
%Remove options already at OPTI default
options.ipopt = optiRemoveDefaults(options.ipopt,ipoptset());
%Add OPTISET options
options.ipopt.print_level = dispLevel(opts.display);
options.ipopt.max_iter = opts.maxiter;
options.ipopt.max_cpu_time = opts.maxtime;
options.ipopt.tol = opts.tolrfun;        

% Get x0 (or try and make one)
if(isfield(prob,'x0') && ~isempty(prob.x0))
    x0 = prob.x0;
    nlprob.x0 = prob.x0; %save it
elseif(prob.sizes.ndec)
    x0 = opti_guessX0(prob,prob.sizes.ndec);
elseif(~isempty(prob.A))
    x0 = opti_guessX0(prob,size(prob.A,2));
elseif(~isempty(prob.nljacstr))
    ndec = size(prob.nljacstr(),2);
    x0 = opti_guessX0(prob,ndec);
elseif(~isempty(prob.Hstr))
    ndec = size(prob.Hstr(),2);
    x0 = opti_guessX0(prob,ndec);
else
    str = sprintf(['OPTI cannot determine the number of decision variables based on the arguments supplied.\n'...
           'Based on this, error checking cannot be performed, and the function cannot continue.\n\nPlease provide x0 to optiprob / opti in order to solve this problem.']);
    error(str); %#ok<SPERR> it doesn't take \n!
end       

%Check obj + grad functions are ok (will segment error otherwise!)
if(isempty(prob.fun))
    error('Your nonlinear objective function is empty - either you are trying to solve a linear problem with a nonlinear solver, or you have forgotten to add your objective (both are errors!)');
end
try        
    f = prob.fun(x0);  %#ok<NASGU>
catch ME
    fprintf(2,'There was an error when running a test objective function call. Please ensure this function exists and runs without error.\n\n');
    rethrow(ME);
end
if(isempty(prob.f))
    error('The gradient function is empty. Please rebuild the OPTI object specifying IPOPT as the solver');
end
try
    g = prob.f(x0); %#ok<NASGU>
catch ME
    fprintf(2,'There was an error when running a test objective gradient function call. Please ensure this function exists and runs without error.\n\n');
    rethrow(ME);
end

%Objective Callback Functions + Hessian
funcs.objective = prob.fun;
funcs.gradient = prob.f;
if(isempty(prob.H))    
    options.ipopt.hessian_approximation = 'limited-memory';
else  
    %Generate a test Hessian
    testH = optiGenTestHessian(prob,x0);    
    %Check we have a tril matrix
    if(any(any(triu(testH,1) ~= 0))); tl = 1; else tl = 0; end
    %Check it is sparse
    if(~issparse(testH)); sp = 1; else sp = 0; end
    %Correct Hessian based on what we found
    switch(nargin(prob.H))
        case 1
            if(tl && sp)            
                funcs.hessian = @(x,sigma,lambda) sigma.*sparse(tril(prob.H(x)));
            elseif(tl)
                funcs.hessian = @(x,sigma,lambda) sigma.*tril(prob.H(x));        
            elseif(sp)
                funcs.hessian = @(x,sigma,lambda) sigma.*sparse(prob.H(x));
            else
                funcs.hessian = @(x,sigma,lambda) sigma.*prob.H(x);
            end
        case 2
            if(tl && sp)            
                funcs.hessian = @(x,sigma,lambda) sigma.*sparse(tril(prob.H(x,lambda)));
            elseif(tl)
                funcs.hessian = @(x,sigma,lambda) sigma.*tril(prob.H(x,lambda));        
            elseif(sp)
                funcs.hessian = @(x,sigma,lambda) sigma.*sparse(prob.H(x,lambda));
            else
                funcs.hessian = @(x,sigma,lambda) sigma.*prob.H(x,lambda);
            end
        case 3
            if(tl && sp)            
                funcs.hessian = @(x,sigma,lambda) sparse(tril(prob.H(x,sigma,lambda)));
            elseif(tl)
                funcs.hessian = @(x,sigma,lambda) tril(prob.H(x,sigma,lambda));        
            elseif(sp)
                funcs.hessian = @(x,sigma,lambda) sparse(prob.H(x,sigma,lambda));
            else
                funcs.hessian = prob.H;
            end
    end    
    %Show warnings
    if(tl && (warn > 1))
        optiwarn('opti_ipopt:hess','Conv IPOPT - The Hessian should be Symmetric TRIL, forcing to TRIL');
    end    
    if(sp && (warn > 1))
        optiwarn('opti_ipopt:hss','Conv IPOPT - The Hessian should return a sparse matrix - correcting');
    end 
    
    %Hessian Structure must be supplied
    if(isempty(prob.Hstr))
        if(warn > 1)
            optiwarn('opti_ipopt:hess','Conv IPOPT - You have not supplied the Hessian Structure of the Objective, Assuming tril(ones(n,n))');
        end
        %Assume all ones
        Hstr = sparse(tril(ones(length(x0))));
        funcs.hessianstructure = @() Hstr;
    else
        try
            testHstr = prob.Hstr();
        catch ME
            fprintf(2,'There was an error when running a test Hessian Structure function call. Please ensure this function exists and runs without error.\n\n');
            rethrow(ME);
        end
        %Check we have a tril matrix
        if(any(any(triu(testHstr,1) ~= 0))); tl = 1; else tl = 0; end        
        if(~issparse(testHstr)); sp = 1; else sp = 0; end
       
        if(tl && sp)            
            funcs.hessianstructure = @() sparse(tril(prob.Hstr()));
        elseif(tl)
            funcs.hessianstructure = @() tril(prob.Hstr());        
        elseif(sp)
            funcs.hessianstructure = @() sparse(prob.Hstr());
        else
            funcs.hessianstructure = prob.Hstr;
        end

        if(sp && (warn > 1))
            optiwarn('opti_ipopt:sp','Conv IPOPT - The Objective Hessian Structure should return a sparse matrix - correcting');
        end
    end
end

%Constraint Callback Functions
if(~isempty(prob.nlcon))
    %Test constraint function
    try
        g = prob.nlcon(x0); %#ok<NASGU>
    catch ME
        fprintf(2,'There was an error when running a test Constraint function call. Please ensure this function exists and runs without error.\n\n');
        rethrow(ME);
    end    
    funcs.constraints = prob.nlcon;
    %Get sample jacobian back
    try
        testJ = prob.nljac(x0);
    catch ME
        fprintf(2,'There was an error when running a test Jacobian function call. Please ensure this function exists and runs without error.\n\n');
        rethrow(ME);
    end
    %transpose decision variable vector (seems odd..?)
    if(size(testJ,2) ~= length(x0)); xt = 1; else xt = 0; end
    %ipopt requires a sparse jacobian
    if(~issparse(testJ)); sp = 1; else sp = 0; end
    %Correct jacobian based on what we found
    if(xt && sp)            
        funcs.jacobian = @(x) sparse(prob.nljac(x'));
    elseif(xt)
        funcs.jacobian = @(x) prob.nljac(x'); 
    elseif(sp)
        funcs.jacobian = @(x) sparse(prob.nljac(x));
    else
        funcs.jacobian = prob.nljac;
    end
    if(sp && (warn > 1))
        optiwarn('opti_ipopt:sp','Conv IPOPT - The Constraint Jacobian should return a sparse matrix - correcting');
    end        
    if(isempty(prob.nljacstr))
        if(warn > 1)
            optiwarn('opti_ipopt:nlcon','Conv IPOPT - You have not supplied the Jacobian Structure of the Constraints, Assuming ones(m,n)');
        end
        %Assume all ones
        Jstr = sparse(ones(length(prob.cl),length(x0)));
        funcs.jacobianstructure = @() Jstr;
    else
        %Get sample jacobian structure back
        try
            testJstr = prob.nljacstr();
        catch ME
            fprintf(2,'There was an error when running a test Jacobian structure function call. Please ensure this function exists and runs without error.\n\n');
            rethrow(ME);
        end
        if(~issparse(testJstr))
            if(warn > 1)
                optiwarn('opti_ipopt:sp','Conv IPOPT - The Constraint Jacobian Structure should return a sparse matrix - correcting');
            end
            funcs.jacobianstructure = @() sparse(prob.nljacstr());
        else
            funcs.jacobianstructure = prob.nljacstr;
        end            
    end
    
end

%Iteration Callback Function
if(isfield(opts,'iterfun') && ~isempty(opts.iterfun))
    funcs.iterfunc = @(i,f,v) iCall(i,f,v,opts.iterfun);
end

%Linear Constraints
if(~isempty(prob.A))
    if(~issparse(prob.A))
        if(warn > 1)
            optiwarn('opti_ipopt:sp','Conv IPOPT - The Linear Constraint matrix A should be a sparse matrix - correcting');
        end
        prob.A = sparse(prob.A);
    end
    options.A = prob.A;
    options.rl = prob.rl;
    options.ru = prob.ru;
end

%Save functions structure
nlprob.funcs = funcs;
%Save Options structure
nlprob.options = options;    

function  print_level = dispLevel(lev)
%Return IPOPT compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 5;
    case 'final'
        print_level = 3;
end

function stop = iCall(iter,fval,vars,iterfun)
%Convert OPTI callback to IPOPT callback
stop = ~iterfun(iter,fval,vars);


