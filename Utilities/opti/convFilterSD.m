function mprob = convFilterSD(prob,opts)
%CONVFILTERSD Convert OPTI problem to FilterSD Nonlinear problem
%
%   mprob = convFilterSD(prob,opts)

%   Copyright (C) 2013 Jonathan Currie (IPL)

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Make sure we have an NLP
if(~any(strcmpi(prob.type,{'NLP','UNO'})))
    error('You can only convert UNOs and NLPs to FILTERSD format!');
end

%Get Warning Level
if(strcmpi(opts.warnings,'all'))
    warn = 2;
elseif(strcmpi(opts.warnings,'critical'))
    warn = 1;
else
    warn = 0;
end

%Convert Row Based Linear Constraints to Row Based Nonlinear Constraints
prob = rowlin2nl(prob,1,warn);
  
%Configure Problem
mprob.lb = prob.lb;
mprob.ub = prob.ub;
mprob.cl = prob.cl;
mprob.cu = prob.cu;

%Add OPTISET options
mprob.options.display = opts.display;
mprob.options.maxiter = opts.maxiter;
mprob.options.maxfeval = opts.maxfeval;
mprob.options.maxtime = opts.maxtime;
mprob.options.rgtol = opts.tolafun;
mprob.options.iterfun = opts.iterfun;

% Get x0 (or try and make one)
if(isfield(prob,'x0') && ~isempty(prob.x0))
    x0 = prob.x0;
    mprob.x0 = prob.x0; %save it
elseif(prob.sizes.ndec)
    x0 = opti_guessX0(prob,prob.sizes.ndec);
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

%Check obj + grad functions are ok (may segment error otherwise!)
try
    f = prob.fun(x0);  %#ok<NASGU>
catch ME
    fprintf(2,'There was an error when running a test objective function call. Please ensure this function exists and runs without error.\n\n');
    rethrow(ME);
end
if(isempty(prob.f))
    error('The gradient function is empty. Please rebuild the OPTI object specifying FilterSD as the solver');
end
try
    g = prob.f(x0); %#ok<NASGU>
catch ME
    fprintf(2,'There was an error when running a test gradient function call. Please ensure this function exists and runs without error.\n\n');
    rethrow(ME);
end

%Objective Callback Functions
mprob.fun = prob.fun;
mprob.grad = prob.f;

%Constraint Callback Functions
if(~isempty(prob.nlcon))
    %Test constraint function
    try
        g = prob.nlcon(x0); %#ok<NASGU>
    catch ME
        fprintf(2,'There was an error when running a test Constraint function call. Please ensure this function exists and runs without error.\n\n');
        rethrow(ME);
    end    
    %Get sample jacobian back
    try
        testJ = prob.nljac(x0);
    catch ME
        fprintf(2,'There was an error when running a test Jacobian function call. Please ensure this function exists and runs without error.\n\n');
        rethrow(ME);
    end
    %Check for sparse Jacobian
    if(issparse(testJ))
        %Check for Jac Structure
        if(isempty(prob.nljacstr))
            if(warn > 1)
                optiwarn('opti_filtersd:jacstr','Conv FILTERSD - You have not supplied the Jacobian Structure of the Constraints, Assuming ones(m,n)');
            end
            %Assume all ones
            Jstr = sparse(ones(length(prob.cl),length(x0)));
            mprob.nljacstr = @() Jstr;
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
                    optiwarn('opti_filtersd:sp','Conv FILTERSD - The Constraint Jacobian Structure should return a sparse matrix - correcting');
                end
                mprob.nljacstr = @() sparse(prob.nljacstr());
            else
                mprob.nljacstr = prob.nljacstr;
            end            
        end
    else
        mprob.nljacstr = []; %no structure as is dense
    end
    
    %Constraint callback functions
    mprob.nlcon = prob.nlcon;
    mprob.nljac = prob.nljac;
else
    mprob.nlcon = [];
    mprob.nljac = [];
    mprob.nljacstr = [];
end


