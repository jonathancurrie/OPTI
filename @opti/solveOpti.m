function [x,fval,exitflag,info] = solveOpti(optObj,x0)
%SOLVE an OPTI object
%
%   Called By opti Solve

%   Copyright (C) 2011-2013 Jonathan Currie (IPL)

%Allocate input args
prob = optObj.prob; 
opts = optObj.opts;
nl = optObj.nlprob;

%Check for solving via AMPL alone
if(prob.ampl.useASL && ~isempty(prob.ampl.path))
    switch(opts.solver)
        case 'scip'
            [x,fval,exitflag,info] = opti_scipasl(prob.ampl.path,opts);
    end
    return
end

%Check initial point
if(~isempty(x0))
    if(numel(x0) ~= optObj.prob.sizes.ndec), error('x0 is not the correct length! Expected %d x 1',optObj.prob.sizes.ndec);end
    if(size(x0,2) > 1), x0 = x0'; end     
    if(any(isnan(x0))), error('Initial Guess (x0) Contains NaN!'); end
    if(any(isinf(x0))), error('Initial Guess (x0) Contains Inf!'); end
    %Ensure dense starting guess
    x0 = full(x0);
    %Save x0 to all problem structures
    if(~isempty(prob)), prob.x0 = x0; end
    if(~isempty(nl)), nl.x0 = x0; end
else
    if(~isempty(nl) && isfield(nl,'x0') && ~isempty(nl.x0))
        if(any(isnan(nl.x0))), error('Initial Guess (x0) Contains NaN!'); end
        if(any(isinf(nl.x0))), error('Initial Guess (x0) Contains Inf!'); end
    elseif(~isempty(prob) && ~isempty(prob.x0))
        if(any(isnan(prob.x0))), error('Initial Guess (x0) Contains NaN!'); end
        if(any(isinf(prob.x0))), error('Initial Guess (x0) Contains Inf!'); end
    end
end


%Switch based on problem type
switch(prob.type)
    case 'SLE'
        [x,fval,exitflag,info] = solveSLE(prob,opts);
    case 'LP'
        [x,fval,exitflag,info] = solveLP(prob,nl,opts);
    case 'MILP'
        [x,fval,exitflag,info] = solveMILP(prob,nl,opts);
    case 'BILP'
        [x,fval,exitflag,info] = solveBILP(prob,nl,opts);
    case 'QP'
        [x,fval,exitflag,info] = solveQP(prob,nl,opts);
    case 'QCQP'
        [x,fval,exitflag,info] = solveQCQP(prob,nl,opts);
    case 'MIQP'
        [x,fval,exitflag,info] = solveMIQP(prob,nl,opts);
    case 'MIQCQP'
        [x,fval,exitflag,info] = solveMIQCQP(prob,nl,opts);
    case 'SDP'
        [x,fval,exitflag,info] = solveSDP(prob,opts);
    case 'MISDP'
        [x,fval,exitflag,info] = solveMISDP(prob,opts);
    case 'SNLE'
        [x,fval,exitflag,info] = solveSNLE(nl,opts,optObj.prob.misc.forceSNLE);
    case 'SCNLE'
        [x,fval,exitflag,info] = solveSCNLE(nl,opts,optObj.prob.misc.forceSNLE);
    case {'DNLS','NLS'}
        [x,fval,exitflag,info] = solveNLS(nl,opts);
    case 'UNO'
        [x,fval,exitflag,info] = solveUNO(nl,opts);    
    case 'NLP'
        [x,fval,exitflag,info] = solveNLP(nl,opts);
    case 'MINLP'
        [x,fval,exitflag,info] = solveMINLP(nl,opts);        
    otherwise
        error('Problem Type ''%s'' is not implemented yet',prob.type);
end



function [x,fval,exitflag,info] = solveSLE(prob,opts)
%Solve a System of Linear Equations
fval = [];
exitflag = [];

switch(opts.solver)
    case 'matlab'
        t = tic;
        x = prob.A\prob.b;
        info = matlabInfo([],[],toc(t),'mldivide');
    case 'mumps'
        [x,info] = opti_mumps(prob.A,prob.b);
    otherwise
        error('The solver %s cannot be used to solve a LSE',opts.solver);
end

function [x,fval,exitflag,info] = solveLP(p,nl,opts)
%Solve a Linear Program using a selected solver

switch(opts.solver)       
    case 'cplex'
        [x,fval,exitflag,info] = opti_cplex([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],[],p.x0,opts.solverOpts); 
    case 'csdp'
        [x,fval,exitflag,info] = opti_csdp(p.f,p.A,p.b,p.lb,p.ub,[],p.x0,opts); 
    case 'dsdp'
        [x,fval,exitflag,info] = opti_dsdp(p.f,p.A,p.b,p.lb,p.ub,[],p.x0,opts); 
    case 'mosek'
        [x,fval,exitflag,info] = moseklp(p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.x0,opts.solverOpts);
    case 'glpk'
        [x,fval,exitflag,info] = opti_glpk(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,opts);
    case 'clp'
        [x,fval,exitflag,info] = opti_clp([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,opts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],[],opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'qsopt'
        [x,fval,exitflag,info] = opti_qsopt(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,opts);
    case 'ooqp'
        [x,fval,exitflag,info] = opti_ooqp([],p.f,p.A,p.rl,p.ru,p.Aeq,p.beq,p.lb,p.ub,opts);
    case 'lp_solve'
        [x,fval,exitflag,info] = opti_lpsolve(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,[],opts);
    case 'lipsol'
        [x,fval,exitflag,info] = opti_lipsol(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,opts);
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
    case 'matlab'
        t = tic;
        [x,fval,exitflag,output,lambda] = linprog(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.x0,opts.solverOpts);
        info = matlabInfo(output,lambda,toc(t),'LINPROG');         
    case 'sedumi'
        [x,fval,exitflag,info] = opti_sedumi(p.sdcone.At,p.sdcone.b,p.sdcone.c,p.sdcone.K,opts);        
    otherwise
        error('The Solver %s cannot be used to solve a LP',opts.solver);        
end

function [x,fval,exitflag,info] = solveMILP(p,nl,opts)
%Solve a Mixed Integer Linear Program using a selected solver

switch(opts.solver)
    case 'cplex'
        [x,fval,exitflag,info] = opti_cplex([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],p.x0,opts.solverOpts);    
    case 'mosek'
        [x,fval,exitflag,info] = mosekmilp(p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.x0,opts.solverOpts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'cbc'
        [x,fval,exitflag,info] = opti_cbc([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,p.x0,opts);
    case 'glpk'
        [x,fval,exitflag,info] = opti_glpk(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,opts);
    case 'lp_solve'
        [x,fval,exitflag,info] = opti_lpsolve(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,p.sos,opts);        
    case 'matlab'
        t = tic;
        [x,fval,exitflag,output] = intlinprog(p.f,find(p.int.str ~= 'C'),p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,opts.solverOpts);
        info = matlabInfo(output,[],toc(t),'INTLINPROG'); 
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);    
    otherwise
        error('The Solver %s cannot be used to solve a MILP',opts.solver);
end

function [x,fval,exitflag,info] = solveBILP(p,nl,opts)
%Solve a Linear Program using a selected solver

switch(opts.solver)       
    case 'cplex'
        [x,fval,exitflag,info] = opti_cplex([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],p.x0,opts.solverOpts);
    case 'mosek'
        [x,fval,exitflag,info] = mosekbilp(p.f,p.A,p.rl,p.ru,p.x0,opts.solverOpts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'glpk'
        [x,fval,exitflag,info] = opti_glpk(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,opts);
    case 'lp_solve'
        [x,fval,exitflag,info] = opti_lpsolve(p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.int.str,p.sos,opts);
    case 'cbc'
        [x,fval,exitflag,info] = opti_cbc([],p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,p.x0,opts);
    case 'matlab'
        if(verLessThan('matlab','8.4')) %2014b
            t = tic;
            [x,fval,exitflag,output] = bintprog(p.f,p.A,p.b,p.Aeq,p.beq,p.x0,opts.solverOpts);
            info = matlabInfo(output,[],toc(t),'BINTPROG');  
        else
            t = tic;
            [x,fval,exitflag,output] = intlinprog(p.f,find(p.int.str ~= 'C'),p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,opts.solverOpts);
            info = matlabInfo(output,[],toc(t),'INTLINPROG');
        end
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);    
    otherwise
        error('The Solver %s cannot be used to solve a BILP',opts.solver);        
end


function [x,fval,exitflag,info] = solveQP(p,nl,opts)
%Solve a Quadratic Program using a selected solver

%Check for unconstrained QP
if(~p.iscon)
    %Solve an Unconstrained QP
    t = tic;
    x = -p.H\p.f;
    info = struct('Iterations',[],'Time',toc(t),'Algorithm','MATLAB: mldivide','StatusString',[]);
    fval = 0.5*x'*p.H*x + p.f'*x;
    exitflag = 1;
else
    switch(opts.solver)       
        case 'cplex'
            [x,fval,exitflag,info] = opti_cplex(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],[],p.x0,opts.solverOpts); 
        case 'mosek'
            [x,fval,exitflag,info] = mosekqp(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.x0,opts.solverOpts);
        case 'ooqp'
            [x,fval,exitflag,info] = opti_ooqp(p.H,p.f,p.A,p.rl,p.ru,p.Aeq,p.beq,p.lb,p.ub,opts);
        case 'clp'
            [x,fval,exitflag,info] = opti_clp(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,opts);
        case 'scip'
            [x,fval,exitflag,info] = opti_scip(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],[],opts);
        case 'baron'
            [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
        case 'ipopt'
            [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
        case 'matlab'
            t = tic;
            [x,fval,exitflag,output,lambda,] = quadprog(p.H,p.f,p.A,p.b,p.Aeq,p.beq,p.lb,p.ub,p.x0,opts.solverOpts);
            info = matlabInfo(output,lambda,toc(t),'QUADPROG'); 
        otherwise
            error('The Solver %s cannot be used to solve a QP',opts.solver);        
    end
end

function [x,fval,exitflag,info] = solveQCQP(p,nl,opts)
%Solve a Quadratically Constrained Program using a selected solver

%Form QC structure
qc.Q = p.Q; qc.l = p.l; qc.qrl = p.qrl; qc.qru = p.qru;

switch(opts.solver)      
    case 'cplex'        
        [x,fval,exitflag,info] = opti_cplex(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],qc,p.x0,opts.solverOpts);
    case 'mosek'
        [x,fval,exitflag,info] = mosekqcqp(p.H,p.f,p.A,p.rl,p.ru,p.Q,p.l,p.qrl,p.qru,p.lb,p.ub,p.x0,opts.solverOpts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,[],[],qc,opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
    otherwise
        error('The Solver %s cannot be used to solve a QCQP',opts.solver);        
end

function [x,fval,exitflag,info] = solveMIQP(p,nl,opts)
%Solve a Mixed Integer Quadratic Program using a selected solver

switch(opts.solver)      
    case 'cplex'
        [x,fval,exitflag,info] = opti_cplex(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],p.x0,opts.solverOpts); 
    case 'mosek'
        [x,fval,exitflag,info] = mosekmiqp(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.x0,opts.solverOpts);
    case 'cbc'
        [x,fval,exitflag,info] = opti_cbc(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,p.x0,opts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,[],opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);        
    otherwise
        error('The Solver %s cannot be used to solve a MIQP',opts.solver);        
end

function [x,fval,exitflag,info] = solveMIQCQP(p,nl,opts)
%Solve a Mixed Integer Quadratically Constrained Quadratic Program using a selected solver

%Form QC structure
qc.Q = p.Q; qc.l = p.l; qc.qrl = p.qrl; qc.qru = p.qru;

switch(opts.solver)      
    case 'cplex'
        [x,fval,exitflag,info] = opti_cplex(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,qc,p.x0,opts.solverOpts);
    case 'mosek'
        [x,fval,exitflag,info] = mosekmiqcqp(p.H,p.f,p.A,p.rl,p.ru,p.Q,p.l,p.qrl,p.qru,p.lb,p.ub,p.int.str,p.x0,opts.solverOpts);
    case 'scip'
        [x,fval,exitflag,info] = opti_scip(p.H,p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.int.str,p.sos,qc,opts);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);        
    otherwise
        error('The Solver %s cannot be used to solve a MIQCQP',opts.solver);        
end

function [x,fval,exitflag,info] = solveSDP(p,opts)
%Solve a Semidefinite Programming problem using a selector solver

switch(opts.solver) 
    case 'csdp'
        [x,fval,exitflag,info] = opti_csdp(p.f,p.A,p.b,p.lb,p.ub,p.sdcone,p.x0,opts);
    case 'dsdp'
        [x,fval,exitflag,info] = opti_dsdp(p.f,p.A,p.b,p.lb,p.ub,p.sdcone,p.x0,opts);
    case 'sedumi'
        [x,fval,exitflag,info] = opti_sedumi(p.sdcone.At,p.sdcone.b,p.sdcone.c,p.sdcone.K,opts);  
    case 'mosek'
        [x,fval,exitflag,info] = moseksdp(p.f,p.A,p.rl,p.ru,p.lb,p.ub,p.sdcone,p.x0,opts.solverOpts);
    otherwise
        error('The Solver %s cannot be used to solve a SDP',opts.solver);        
end

function [x,fval,exitflag,info] = solveMISDP(p,opts)
%Solve a Mixed Integer Semidefinite Programming problem using a selector solver
switch(opts.solver) 
    case 'scip'
        [x,fval,exitflag,info] = opti_scipsdp(p.f,p.A,p.b,p.lb,p.ub,p.sdcone,p.int.str,opts);
    otherwise
        error('The Solver %s cannot be used to solve a MISDP',opts.solver);  
end


function [x,fval,exitflag,info] = solveSNLE(nl,opts,force)
%Solve a Nonlinear Equation Problem using a selected solver

%Allows for runtime generation of eq rhs
if(~isfield(nl,'ydata') || isempty(nl.ydata)) 
    nl.ydata = zeros(size(nl.x0));
end

switch(opts.solver)
    case 'matlab'
        t = tic;
        [x,fval,exitflag,output] = fsolve(nl);
        info = matlabInfo(output,[],toc(t),'FSOLVE');
    case 'hybrj'
        [x,fval,exitflag,info] = opti_hybrj(nl.fun,nl.grad,nl.x0,nl.options);
    case 'nl2sol'        
        [x,fval,exitflag,info] = opti_nl2sol(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.options);
    case 'lmder'
        [x,fval,exitflag,info] = opti_lmder(nl.fun,nl.grad,nl.x0,nl.ydata,nl.options);
    case 'mkltrnls'
        [x,fval,exitflag,info] = opti_mkltrnls(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,opts);
    case 'levmar'
        [x,fval,exitflag,info] = opti_levmar(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.A,nl.b,nl.Aeq,nl.beq,nl.options);
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0); 
        if(force), fval = sqrt(sum(nl.funcs.constraints(x).^2)); end
    case 'nomad'
        [x,~,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'nlopt'
        [x,fval,exitflag,info] = opti_nlopt(nl,nl.x0); 
        if(force), fval = sqrt(sum(nl.nlcon(x).^2)); end
    otherwise
        error('The Solver %s cannot be used to solve a SNLE',opts.solver);
end

function [x,fval,exitflag,info] = solveSCNLE(nl,opts,force)
%Solve a Constrained Nonlinear Equation Problem using a selected solver

%Allows for runtime generation of eq rhs
if(~isfield(nl,'ydata') || isempty(nl.ydata)) 
    nl.ydata = zeros(size(nl.x0));
end

switch(opts.solver)
    case 'nl2sol'        
        [x,fval,exitflag,info] = opti_nl2sol(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.options);
    case 'mkltrnls'
        [x,fval,exitflag,info] = opti_mkltrnls(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,opts);
    case 'levmar'
        [x,fval,exitflag,info] = opti_levmar(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.A,nl.b,nl.Aeq,nl.beq,nl.options);
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
        if(force), fval = sqrt(sum(nl.funcs.constraints(x).^2)); end
    case 'nomad'
        [x,fval,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'nlopt'
        [x,fval,exitflag,info] = opti_nlopt(nl,nl.x0); 
        if(force), fval = sqrt(sum(nl.funcs.constraints(x).^2)); end
    case 'pswarm'
        [x,fval,exitflag,info] = opti_pswarm(nl.fun,nl.lb,nl.ub,nl.x0,nl.A,nl.b,nl.options);
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);
        if(force), fval = sqrt(sum(nl.funcs.constraints(x).^2)); end
    otherwise
        error('The Solver %s cannot be used to solve a SCNLE',opts.solver);
end

function [x,fval,exitflag,info] = solveNLS(nl,opts)
%Solve a Nonlinear Least Squares Problem using a selected solver

switch(opts.solver)
    case 'matlab'
        t = tic;
        [x,fval,~,exitflag,output,lambda] = lsqnonlin(nl);
        info = matlabInfo(output,lambda,toc(t),'LSQNONLIN');
    case 'lmder'
        [x,fval,exitflag,info] = opti_lmder(nl.fun,nl.grad,nl.x0,nl.ydata,nl.options);
    case 'levmar'
        [x,fval,exitflag,info] = opti_levmar(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.A,nl.b,nl.Aeq,nl.beq,nl.options);
    case 'mkltrnls'
        [x,fval,exitflag,info] = opti_mkltrnls(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,opts);
	case 'nl2sol'
		[x,fval,exitflag,info] = opti_nl2sol(nl.fun,nl.grad,nl.x0,nl.ydata,nl.lb,nl.ub,nl.options); 
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
    case 'nomad'
        [x,fval,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'nlopt'
        [x,fval,exitflag,info] = opti_nlopt(nl,nl.x0); 
    case 'pswarm'
        [x,fval,exitflag,info] = opti_pswarm(nl.fun,nl.lb,nl.ub,nl.x0,nl.A,nl.b,nl.options);
    otherwise
        error('The Solver %s cannot be used to solve a NLS',opts.solver);
end
%Close CVODEs if opened
if(isfield(opts,'dynamicOpts') && isfield(opts.dynamicOpts,'integrator') && strcmpi(opts.dynamicOpts.integrator,'cvodes'))
    CVodeFree;
end

function [x,fval,exitflag,info] = solveUNO(nl,opts)
%Solve an Unconstrained Nonlinear Problem

switch(opts.solver)
    case 'matlab'
        t = tic;
        [x,fval,exitflag,output] = fminsearch(nl.objective,nl.x0);
        info = matlabInfo(output,[],toc(t),'FMINSEARCH - Simplex');
    case 'filtersd'
        [x,fval,exitflag,info] = opti_filtersd(nl.fun,nl.grad,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nljac,nl.nljacstr,nl.cl,nl.cu,nl.options);
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0);
    case 'nlopt'
        [x,fval,exitflag,info] = opti_nlopt(nl,nl.x0);
    case 'nomad'
        [x,fval,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'm1qn3'
        [x,fval,exitflag,info] = opti_m1qn3(nl.fun,nl.grad,nl.x0,nl.options);
    case 'scip'
        [x,fval,exitflag,info] = opti_scipnl(nl.fun,[],[],[],[],[],[],[],[],[],nl.x0,nl.options);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'gmatlab'
        t = tic;
        [x,fval,exitflag,output] = patternsearch(nl);
        info = gmatlabInfo(output,toc(t),'PATTERNSERACH - Global Direct Search');
    otherwise
        error('The Solver %s cannot be used to solve this type of problem',opts.solver);
end

function [x,fval,exitflag,info] = solveNLP(nl,opts)
%Solve a Nonlinear Program using a selected solver

if(~isfield(nl,'x0')), error('You must supply x0 to solve an NLP'); end

switch(opts.solver)
    case 'matlab'
        t = tic;
        [x,fval,exitflag,output,lambda] = fmincon(nl);
        info = matlabInfo(output,lambda,toc(t),'FMINCON');
    case 'ipopt'
        [x,fval,exitflag,info] = opti_ipopt(nl,nl.x0); 
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);
    case 'filtersd'
        [x,fval,exitflag,info] = opti_filtersd(nl.fun,nl.grad,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nljac,nl.nljacstr,nl.cl,nl.cu,nl.options);
    case 'nlopt'
        [x,fval,exitflag,info] = opti_nlopt(nl,nl.x0); 
    case 'lbfgsb'
        [x,fval,exitflag,info] = opti_lbfgsb(nl.fun,nl.grad,nl.lb,nl.ub,nl.x0,nl.options);
    case 'nomad'
        [x,fval,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'pswarm'
        [x,fval,exitflag,info] = opti_pswarm(nl.fun,nl.lb,nl.ub,nl.x0,nl.A,nl.b,nl.options);
    case 'scip'
        [x,fval,exitflag,info] = opti_scipnl(nl.fun,nl.A,nl.rl,nl.ru,nl.lb,nl.ub,nl.nlcon,nl.cl,nl.cu,nl.int.str,nl.x0,nl.options);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'gmatlab'
        t = tic;
        [x,fval,exitflag,output] = patternsearch(nl);
        info = gmatlabInfo(output,toc(t),'PATTERNSERACH - Global Direct Search');
    otherwise
        error('The Solver %s cannot be used to solve a NLP',opts.solver);
end

function [x,fval,exitflag,info] = solveMINLP(nl,opts)
%Solve a Mixed Integer Nonlinear Program using a selected solver

%Should not happen, but worthwhile to check
if(isempty(nl))
    error('The Solver %s cannot be used to solve a MINLP',opts.solver);
end

if(~isfield(nl,'x0')), error('You must supply x0 to solve an MINLP'); end

switch(opts.solver)
    case 'bonmin'
        [x,fval,exitflag,info] = opti_bonmin(nl,nl.x0);
    case 'scip'
        [x,fval,exitflag,info] = opti_scipnl(nl.fun,nl.A,nl.rl,nl.ru,nl.lb,nl.ub,nl.nlcon,nl.cl,nl.cu,nl.int.str,nl.x0,nl.options);
    case 'baron'
        [x,fval,exitflag,info] = opti_baron(nl,nl.x0);
    case 'nomad'
        [x,fval,exitflag,info] = opti_nomad(nl.fun,nl.x0,nl.lb,nl.ub,nl.nlcon,nl.nlrhs,nl.xtype,nl.options);
    case 'gmatlab'
        t = tic;
        [x,fval,exitflag,output] = ga(nl);
        info = gmatlabInfo(output,toc(t),'GA - Genetic Algorithm');
    otherwise
        error('The Solver %s cannot be used to solve a MINLP',opts.solver);
end


function info = matlabInfo(output,lambda,t,alg)
%Convert Matlab Output to opti info structure

if(isempty(output))
    output.iterations = [];
    output.message = [];
end
if(~isfield(output,'message'))
    output.message = 'This is does not appear to be a MATLAB solver - You have overloaded it?';
end
if(~isfield(output,'algorithm'))
    ALG = 'MATLAB';
else
    ALG = ['MATLAB: ' alg ' - ' output.algorithm];
end
if(isfield(output,'iterations'))
    info = struct('Iterations',output.iterations);
else
    info = struct('Nodes',output.numnodes);
end
info.Time = t;
info.Algorithm = ALG;
if(~isempty(output.message))
    msg = regexp(output.message,'\n','split');
    info.Status = msg{1};
else
    info.Status = '';
end

if(~isempty(lambda))
    info.Lambda = lambda;
end


function info = gmatlabInfo(output,t,alg)
%Convert GMatlab Output to opti info structure

if(isempty(output))
    output.iterations = [];
    output.message = [];
end
if(~isfield(output,'message'))
    output.message = 'This is does not appear to be a MATLAB solver  - You have overloaded it?';
end

if(isfield(output,'iterations'))
    info.Iterations = output.iterations;
elseif(isfield(output,'generations'))
    info.Iterations = output.generations;
else
    info.Iterations = [];
end
info.Time = t;
info.Algorithm = ['GMATLAB: ' alg];
info.Status = output.message;                       