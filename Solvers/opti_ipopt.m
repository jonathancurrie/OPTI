function [x,fval,exitflag,info] = opti_ipopt(nlprob,x0)
%OPTI_IPOPT Solve a NLP using IPOPT
%
%   [x,fval,exitflag,info] = opti_ipopt(nlprob,x0) solves the nonlinear
%   program min f(x) subject to linear and nonlinear constraints using
%   IPOPT. nlprob is supplied in nl opti format (i.e. from convIpopt())
%   and x0 is the initial solution guess.
%
%   THIS IS A WRAPPER FOR IPOPT
%   See supplied Eclipse Public License

%   Copyright (C) 2011-2013 Jonathan Currie (Control Engineering)

t = tic;
%Check required fields
if(~isfield(nlprob,'funcs') || ~isfield(nlprob,'options'))
    error('You must use convert(optiObj) or solve(optiObj) to generate a problem structure for this function');
end
%Ensure we have a starting guess
if(nargin < 2 || isempty(x0))
    if(isfield(nlprob,'x0') && ~isempty(nlprob.x0))
        x0 = nlprob.x0;
    else
        error('You must supply x0 to use ipopt!');
    end
end

% Set optiver for version comparison
nlprob.options.optiver = optiver;

% Run IPOPT
[x,output] = ipopt(x0,nlprob.funcs,nlprob.options);

%Collect Results
fval = nlprob.funcs.objective(x);
info.Iterations = output.iter;
info.FuncEvals = output.eval;
info.Time = toc(t);
info.Algorithm = 'IPOPT: Interior Point NL Solver';

switch(output.status)
    case 0
        info.Status = 'Success';
        exitflag = 1;
    case 1
        info.Status = 'Solved to Acceptable Level';
        exitflag = 1;
    case 2
        info.Status = 'Infeasible';
        exitflag = -1;
    case 3
        info.Status = 'Search Direction Becomes Too Small';
        exitflag = -2;
    case 4
        info.Status = 'Diverging Iterates';
        exitflag = -2;
    case 5
        info.Status = 'User Exit';
        exitflag = -5;
    case 6
        info.Status = 'Feasible Point Found';
        exitflag = 1;
        
    case -1
        info.Status = 'Exceeded Iterations';
        exitflag = 0;
    case -2
        info.Status = 'Restoration Failed';
        exitflag = -3;
    case -3
        info.Status = 'Error in Step Computation';
        exitflag = -3;
    case -4 
        info.Status = 'Max Time Exceeded';
        exitflag = 0;
    case -10
        info.Status = 'Not Enough Degrees Of Freedom';
        exitflag = -3;
    case -11
        info.Status = 'Invalid Problem Definition';
        exitflag = -4;
    case -12
        info.Status = 'Invalid Option';
        exitflag = -4;
    case -13
        info.Status = 'Invalid Number Detected';
        exitflag = -3;
    case -102
        info.Status = 'Insufficient Memory';
        exitflag = -3;
    otherwise        
        info.Status = 'IPOPT Error';
        exitflag = -3;
end

%Assign Lambda
lam = output.lambda;
if(~isempty(lam))
    info.Lambda.ineqlin = lam(nlprob.options.ineq);
    info.Lambda.eqlin = lam(nlprob.options.eq);
    info.Lambda.ineqnonlin = lam(nlprob.options.nlineq);
    info.Lambda.eqnonlin = lam(nlprob.options.nleq);
end
info.Lambda.upper = output.zu;
info.Lambda.lower = output.zl;
