function mprob = convNlopt(prob,opts)
%CONVNLOPT Convert OPTI problem to NLOPT Nonlinear problem
%
%   mprob = convNlopt(prob,opts)

%   Copyright (C) 2011 Jonathan Currie (I2C2)

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Get Warning Level
warn = optiWarnLevel(opts.warnings);

%Check for and process linear constraints
prob = genlin2nl(prob,0,warn);

%Read in passed nlopt options
setOpts = nloptset(opts.solverOpts);

% Set NLOPT options 
mprob.options.tolafun    = opts.tolafun;
mprob.options.tolrfun    = opts.tolrfun;
mprob.options.maxtime    = opts.maxtime;
mprob.options.maxfeval   = opts.maxfeval;
mprob.options.iterfun    = opts.iterfun;
mprob.options.display    = dispLevel(opts.display);
mprob.options.algorithm  = nloptSolver(setOpts.algorithm);
%Set Suboptimizer options
mprob.options.local_optimizer.algorithm = nloptSolver(setOpts.subalgorithm);
mprob.options.local_optimizer.tolrfun = setOpts.subtolrfun;    
mprob.options.local_optimizer.tolafun = setOpts.subtolafun;    
mprob.options.local_optimizer.maxtime = setOpts.submaxtime;
mprob.options.local_optimizer.maxfeval = setOpts.submaxfeval;

% Get x0 (or try and make one for testing only)
if(isfield(prob,'x0') && ~isempty(prob.x0))
    x0 = prob.x0;
elseif(prob.sizes.ndec)
    x0 = opti_guessX0(prob,prob.sizes.ndec);
else
    str = sprintf(['OPTI cannot determine the number of decision variables based on the arguments supplied.\n'...
           'Based on this, error checking cannot be performed, and the function cannot continue.\n\nPlease provide x0 to optiprob in order to solve this problem.']);
    error(str); %#ok<SPERR> it doesn't take \n!
end
                                          
%Configure Problem Bounds
mprob.lb = prob.lb;
mprob.ub = prob.ub; 

% Check NLOPT Algorithm Selection
[~,ineq,eq,uncon,~,deriv] = nloptSolver(mprob.options.algorithm);
% Check gradient settings
if(deriv && prob.numdif.grad && warn && mprob.options.algorithm~=36) %derivative based solver but with numerical difference scheme (bad)
    optiwarn('opti_nlopt:numdiff',['Conv NLOPT - The currently selected NLOPT algorithm (%s) requires a gradient function '...
                                  'but you have not supplied one. Try using a derivative free algorithm instead '...
                                  '(check nloptSolver)'],upper(setOpts.algorithm));
end
if(~deriv && ~prob.numdif.grad && ~isempty(prob.f) && warn && mprob.options.algorithm~=36) %not derivative based solver, but supplied gradient (also bad)
    optiwarn('opti_nlopt:numdiff',['Conv NLOPT - The currently selected NLOPT algorithm (%s) does not require a gradient function '...
                                  'but you have supplied one. Try using a derivative enabled algorithm instead '...
                                  '(check nloptSolver)'],upper(setOpts.algorithm));
end
% Check if unconstrained with not an unconstrained solver
if(~prob.iscon && ~uncon)
    if(prob.numdif.grad) %check if we are using numerical derivative information (bad? we will use deriv free opt with nlopt if estimated)
        alg = 'LN_PRAXIS'; %seems to work well - you are welcome to change the default!
    else
        alg = 'LD_LBFGS'; %also seems to work ok
    end    
    if(warn)
        optiwarn('opti_nlopt:uncon','Conv NLOPT - The currently selected NLOPT algorithm (%s) cannot solve unconstrained problems, using %s instead',...
                upper(nloptSolver(mprob.algorithm)),alg);
    end
    [~,ineq,eq,~,~,deriv] = nloptSolver(alg); %recheck
    mprob.algorithm = nloptSolver(alg);
end
%Check if nonlinearly inequal constrained but using solver that can't support them
if(~ineq && prob.sizes.nnlineq)
    if(deriv) %using a derivate based solver above
        alg = 'LD_SLSQP';
    else
        alg = 'LN_AUGLAG';
    end
    if(warn)
        optiwarn('opti_nlopt:uncon','Conv NLOPT - The currently selected NLOPT algorithm (%s) cannot solve problems with nonlinear inequality constraints, using %s instead',...
                upper(nloptSolver(mprob.algorithm)),alg);
    end
    [~,~,eq,~,~,deriv] = nloptSolver(alg); %recheck
    mprob.algorithm = nloptSolver(alg);
end
%Check if nonlinearly equal constrained but using a solver that can't support them
if(~eq && prob.sizes.nnleq)
    if(deriv) %using a derivate based solver above
        alg = 'LD_SLSQP';
    else
        alg = 'LN_AUGLAG';
    end
    if(warn && isfield(mprob,'algorithm'))
        optiwarn('opti_nlopt:uncon','Conv NLOPT - The currently selected NLOPT algorithm (%s) cannot solve problems with nonlinear equality constraints, using %s instead',...
                upper(nloptSolver(mprob.algorithm)),alg);
    end
    mprob.algorithm = nloptSolver(alg);
end

%Final Properties Get & Save
[~,~,~,uncon,glob,deriv,subprob] = nloptSolver(mprob.options.algorithm);
mprob.alg_prop = struct('uncon',uncon,'glob',glob,'deriv',deriv,'subopt',subprob);

%If we aren't using local optimizer, remove the field
if(~subprob && isfield(mprob,'local_optimizer'))
    mprob = rmfield(mprob,'local_optimizer');
end

%Objective & Gradient Callback Function
mprob.objective = prob.fun;
mprob.gradient = prob.f;

%Constraints
mprob.nlcon = prob.nlcon;
mprob.nlrhs = prob.nlrhs;
mprob.nle = prob.nle;

%Check jac is not sparse
if(~isempty(prob.nljac))    
    testJ = prob.nljac(x0);
    if(issparse(testJ))            
        if(warn > 1)
            optiwarn('opti_nlopt:sp','Conv NLOPT - The Constraint Jacobian Function should return a dense matrix - correcting: [full(nljac)]');
        end
        mprob.nljac = @(x) full(prob.nljac(x));
    else
        mprob.nljac = prob.nljac;
    end
end