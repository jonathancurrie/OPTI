function mprob = convBaron(prob,opts)
%CONVBARON Convert OPTI problem to BARON Nonlinear problem
%
%   mprob = convBaron(prob,opts)

%   Copyright (C) 2013 Jonathan Currie (I2C2)

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

%Assign Common Args
mprob.A = prob.A;
mprob.rl = prob.rl;
mprob.ru = prob.ru;
mprob.lb = prob.lb;
mprob.ub = prob.ub;
if(isfield(prob,'int') && isfield(prob.int,'str'))
    mprob.xtype = prob.int.str;
end

%Check if the user wants to use BARON LP/QP/QCQP solver
switch(lower(prob.type))
    case {'lp','milp','bilp','qp','qcqp','miqp','miqcqp'}
        if(warn > 1)
            optiwarn('OPTI:QCQP2NLP','Quadratic Program converted to (MI)NLP to suit BARON interface');
        end
        [prob,opts] = QCQP2NLP(prob,opts);
    case {'uno','nlp','minlp'}
        %nothing
    otherwise
        error('Problem type ''%s'' is not supported by BARON',prob.type);
end

%Assign remaining args
mprob.fun = prob.fun;
mprob.nlcon = prob.nlcon;
mprob.cl = prob.cl;
mprob.cu = prob.cu;
mprob.x0 = prob.x0;

% Set BARON options (optiset options override common ones)
%Build baronset struct, adding any passed options (ignored by baronset if not recognised)
options = baronset(opts.solverOpts);
%Add OPTISET options
options.epsa = opts.tolafun;
options.epsr = opts.tolrfun;
options.inttol = opts.tolint;
options.maxiter = opts.maxiter;
options.maxtime = opts.maxtime;    
options.prlevel = dispLevel(opts.display);
%Save Options Structure
mprob.options = options;    

function  print_level = dispLevel(lev)
%Return IPOPT compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 1;
    case 'final'
        print_level = 1;
end




