function mprob = convBonmin(prob,opts,warn)
%CONVBONMIN Convert OPTI problem to BONMIN Mixed Integer Nonlinear problem
%
%   mprob = convBonmin(prob,opts)

%   Copyright (C) 2011 Jonathan Currie (I2C2)

%COMPATIBLE CPLEX VERSION
CPLEXVER = '12.5.1.0';

%Ensure all args passed
if(nargin < 2)
    error('You must supply both the problem + options');
end
if(~isstruct(prob) || ~isstruct(opts))
    error('Both prob + opts must be structures!');
end

%Ensure we have defaults for bonmin
if(isempty(opts.solverOpts))
    opts.solverOpts = bonminset;
end

%First convert to IPOPT problem (used as the internal relaxed solver)
mprob = convIpopt(prob,opts);
%Remove options from IPOPT structure we don't need
mprob.options.ipopt = rmfield(mprob.options.ipopt,'max_cpu_time'); %don't want this set in the local solver
mprob.options.ipopt = rmfield(mprob.options.ipopt,'print_level'); %don't want this set in the local solver
%Remove options already at OPTI default (appears BONMIN sets some internally different from docs)
mprob.options.ipopt = optiRemoveDefaults(mprob.options.ipopt,ipoptset('bonmin'));

%Set BONMIN options
mprob.options.bonmin = bonminset(opts.solverOpts,'noIpopt'); %this will remove ipopt options from the bonminset structure
%Remove options set elsewhere
mprob.options.bonmin = rmfield(mprob.options.bonmin,'var_lin'); 
mprob.options.bonmin = rmfield(mprob.options.bonmin,'cons_lin'); 
%Remove options already at OPTI default (appears BONMIN sets some internally different from docs)
mprob.options.bonmin = optiRemoveDefaults(mprob.options.bonmin,bonminset());

%If user has specified Cplex, ensure we can load it
if(isfield(mprob.options.bonmin,'milp_solver') && strcmpi(mprob.options.bonmin.milp_solver,'cplex'))
    try
        %Forces to load Cplex on current platform (x86 or x64) if MATLAB interface present
        c = Cplex; %#ok<NASGU>
    catch %#ok<CTCH>
       %nothing to do 
    end
    try
        %The following will fail with a "cannot load module" if it cannot find correct version of Cplex DLL on user's PC
        v = bonminCplex(); %#ok<NASGU>
    catch %#ok<CTCH>
        if(warn)
            optiwarn('opti:bonmin_Cplex','Could not load BONMIN with CPLEX - ensure CPLEX (v%s) is installed and licensed on your PC. Using CBC for now.',CPLEXVER);
        end
        mprob.options.bonmin = rmfield(mprob.options.bonmin,'milp_solver'); 
    end
end

%Convert OPTISET options to equivalent BONMIN ones
mprob.options.bonmin.node_limit = opts.maxnodes;
mprob.options.bonmin.time_limit = opts.maxtime;
mprob.options.bonmin.integer_tolerance = opts.tolint;

%Set Display Level
if(strcmpi(opts.display,'iter'))
    mprob.options.display = 1;
else
    mprob.options.display = 0;
end

%Setup Integer Vars
mprob.options.var_type = prob.int.ind; %uses -1, 0, 1 format
%Setup Nonlinear Vars
if(~isfield(opts.solverOpts,'var_lin') || isempty(opts.solverOpts.var_lin) || any(isnan(opts.solverOpts.var_lin)))
    mprob.options.var_lin = zeros(prob.sizes.ndec,1); %for now assuming all nonlinear
else
    varlin = opts.solverOpts.var_lin;
    if(length(varlin) ~= prob.sizes.ndec)
        error('The decision variable linearity vector is the wrong size! Expected %d x 1',prob.sizes.ndec);
    else
        mprob.options.var_lin = varlin;
    end
end
%Setup Nonlinear Constraints Linearity
if(prob.sizes.nrow) %bug fix when bonmin solves with double-sided linear row constraints
    nlincon = prob.sizes.nrow;
else
    nlincon = prob.sizes.nineq + prob.sizes.neq;
end
if(prob.sizes.nnlrow) %bug fix when bonmin solves with double-sided nonlinear row constraints
    nnlcon = prob.sizes.nnlrow + prob.sizes.nqc;
else
    nnlcon = prob.sizes.nnlineq + prob.sizes.nnleq + prob.sizes.nqc;
end
if(~isfield(opts.solverOpts,'cons_lin') || isempty(opts.solverOpts.cons_lin) || any(isnan(opts.solverOpts.cons_lin)))
   mprob.options.cons_lin = [zeros(nnlcon,1); ones(nlincon,1)]; %augment both kinds
else
    conslin = opts.solverOpts.cons_lin;
    if(length(conslin) ~= (nlincon+nnlcon))
        error('The constraint linearity vector is the wrong size! Expected %d x 1',nlincon+nnlcon);
    else
        mprob.options.cons_lin = conslin;
    end
end