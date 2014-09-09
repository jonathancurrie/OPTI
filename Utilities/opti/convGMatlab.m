function mprob = convGMatlab(prob,opts)
%CONVGMATLAB Convert OPTI problem to MATLAB Global Optimization Toolbox problem
%
%   mprob = convGMatlab(prob,opts)

%   Copyright (C) 2012 Jonathan Currie (I2C2)

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Make sure we have an NLP
if(~any(strcmpi(prob.type,{'NLP','UNO','MINLP'})))
    error('You can only convert UNOs, NLPs and MINLPs to MATLAB Global format!');
end

%Starting Guess
if(isfield(prob,'x0'))
    mprob.x0 = prob.x0;
end

warn = strcmpi(opts.warnings,'all');

%Problem dependent args
switch(upper(prob.type))
    case 'UNO'
        mprob.solver = 'patternsearch';
        mprob.options = psoptimset(opts.solverOpts,'TimeLimit',opts.solverOpts.MaxTime);
        mprob.objective = prob.fun;       
    
    case 'NLP'       
        mprob.solver = 'patternsearch'; 
        mprob.options = psoptimset(opts.solverOpts,'TimeLimit',opts.solverOpts.MaxTime);
        %Linear stuff
        mprob.Aineq = prob.A;
        mprob.bineq = prob.b;
        mprob.Aeq = prob.Aeq;
        mprob.beq = prob.beq;
        mprob.lb = prob.lb;
        mprob.ub = prob.ub;
        %Objective
        mprob.objective = prob.fun;       
        %Check for NL Constraints
        if(~isempty(prob.nlcon))
            %Get Constraint Types
            max_in = prob.nle == 1;
            min_in = prob.nle == -1;
            eq = prob.nle == 0;
            mprob.nonlcon = @(x) nlCon(x,prob.nlcon,prob.nlrhs,max_in,min_in,eq);
        end
    case 'MINLP'
        mprob.solver = 'ga'; 
        mprob.options = psoptimset(opts.solverOpts,'TimeLimit',opts.solverOpts.MaxTime);
        %Linear stuff
        mprob.Aineq = prob.A;
        mprob.bineq = prob.b;
        mprob.lb = prob.lb;
        mprob.ub = prob.ub;
        if(~isempty(prob.beq) && warn)
            warning('opti_gmatlab:ga_eq','The MATLAB GA Integer Solver does not accept linear equalities, ignoring');
        end
        %Objective
        mprob.fitnessfcn = prob.fun;     
        mprob.nvars = prob.sizes.ndec;
        %Check for NL Constraints
        if(~isempty(prob.nlcon))
            %Get Constraint Types
            max_in = prob.nle == 1;
            min_in = prob.nle == -1;
            eq = prob.nle == 0;
            if(~isempty(eq) && warn)
                warning('opti_gmatlab:ga_eq','The MATLAB GA Integer Solver does not accept nonlinear equalities, ignoring');
            end
            mprob.nonlcon = @(x) nlCon(x,prob.nlcon,prob.nlrhs,max_in,min_in,eq);
        end
        %Integer Constraints
        if(any(prob.int.ind))
            mprob.intcon = prob.int.idx;
        end
        
    otherwise
        error('Not implemented yet');
end

%Check for iterfun
if(isfield(opts,'iterfun') && ~isempty(opts.iterfun))
    mprob.options.OutputFcns = @(optimV,o,f) mCall(optimV,o,f,opts.iterfun);
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


function [stop,o,oc] = mCall(optimValues,o,state,fun)
% Handle to convert MATLAB callback to OPTI callback

stop = false;
oc = false;
switch state
    case 'iter'
        stop = fun(optimValues.iteration,optimValues.fval,optimValues.x);
        drawnow;
end
