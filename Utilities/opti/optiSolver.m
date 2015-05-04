function ret = optiSolver(instr,err)
%OPTISOLVER  See if a solver is installed on your PC and details about it
%
%   optiSolver() prints a list of OPTI interfaced solvers, their version
%   numbers and and whether they are available on your PC.
%
%   optiSolver('matrix') prints a solver vs problem type matrix for
%   determining which solvers and setup for which problems.
%
%   optiSolver('config') prints a solver vs config matrix for
%   determining which options are configured for which solvers.
%
%   optiSolver(type) prints a list of solvers for a given problem, 'type',
%   together with a table of constraint and derivative information. See
%   below for heading information. Examples are 'lp', 'milp', etc.
%
%   ok = optiSolver(solver) returns true if the passed solver is available
%   on your PC. Examples are 'clp', 'ipopt', etc.
%
%   optiSolver(solver) displays information about the solver (if
%   available) and configuration specifics (no output args).
%
%   solver = optiSolver(bestSolver) returns the best available solver as a
%   string for a given problem type. Examples are 'best_lp', 'best_milp',
%   etc.
%
%   Problem Type table headings are as follows:
%
%   BD - Bounds
%   LI - Linear Inequality (* indicates a conversion will take place)
%   LE - Linear Equality   (* indicates a conversion will take place)
%   QI - Quadratic Inequality
%   SS - Special Ordered Set
%   NI - Nonlinear Inequality
%   NE - Nonlinear Equality
%   SP - Supports Sparse Matrices
%   D1 - Requires First Derivatives (Gradient) (* indicates optional)
%   D2 - Requires Second Derivatives (Hessian) (* indicates optional)
%   GL - Global Solver

%   Copyright (C) 2011-2015 Jonathan Currie (I2C2)


global SLE LP MILP BILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP DESC PTYPES

%Ordered lists of 'best' solvers
SLE = {'matlab','MUMPS'};
LP = {'CPLEX','clp','mosek','scip','ooqp','glpk','lp_solve','MATLAB','csdp','dsdp','SeDuMi'};
MILP = {'CPLEX','cbc','mosek','scip','glpk','lp_solve','MATLAB'};
BILP = {'CPLEX','cbc','mosek','scip','glpk','lp_solve','MATLAB'};
QP = {'CPLEX','ooqp','clp','mosek','scip','ipopt','MATLAB'}; 
QCQP = {'CPLEX','mosek','scip','ipopt'};
MIQP = {'CPLEX','mosek','scip','bonmin'};
MIQCQP = {'CPLEX','mosek','scip','bonmin'};
SDP = {'csdp','dsdp','SeDuMi'}; %mosek?
MISDP = {}; %'scip'
SNLE = {'nl2sol','hybrj','mkltrnls','lmder','levmar','MATLAB','IPOPT','nomad','nlopt'};
SCNLE = {'nl2sol','mkltrnls','levmar','IPOPT','nomad','pswarm','nlopt'};
NLS = {'nl2sol','mkltrnls','lmder','levmar','MATLAB','IPOPT','nomad','pswarm','nlopt'};
UNO = {'IPOPT','matlab','m1qn3','NLOPT','nomad','scip','GMATLAB'};
NLP = {'IPOPT','MATLAB','FILTERSD','NLOPT','LBFGSB','nomad','PSwarm','baron','scip','GMATLAB'};
MINLP = {'BONMIN','nomad','baron','scip','GMATLAB'};
%Description List
DESC = containers.Map({'cplex','glpk','lp_solve','clp','matlab','mumps','ooqp','ipopt','nlopt','cbc','bonmin','lbfgsb','levmar','lmder','nl2sol','hybrj',...
                       'mkltrnls','pswarm','gmatlab','nomad','scip','m1qn3','filtersd','dsdp','csdp','sedumi','baron','mosek'},...
                      {'IBM ILOG CPLEX','GNU Linear Programming Kit','LP_Solve','COIN-OR Linear Programming',...
                       'Optimization Toolbox: ','a MUltifrontal Massively Parallel sparse direct Solver','Object-Orientated Quadratic Programming',...
                       'Interior-Point Optimization','Nonlinear Optimization','COIN-OR Branch and Cut',...
                       'Basic Open-source Nonlinear Mixed INteger programming','Limited Memory Broyden-Fletcher-Goldfarb-Shanno Bounded Optimization',...
                       'Levenberg-Marquardt Nonlinear Least Squares in C/C++','MINPACK Levenberg-Marquardt','Adaptive Nonlinear Least Squares','MINPACK Powell Hybrid',...
                       'Intel MKL Trust Region Nonlinear Least Squares','Pattern and Particle Swarm Global Optimizer','Global Optimization Toolbox: ',...
                       'Nonlinear Optimization with the MADS algorithm','Solving Constraint Integer Programs','Large-Scale Unconstrained L-BFGS Minimization',...
                       'FilterSD: Filter & Trust Region Nonlinear Optimization','Dual Scaling Interior Point SDP Solver','Predictor-Corrector Primal-Dual SDP Solver',...
                       'Self-Dual-Minimization','Branch And Reduce Optimization Navigator','MOSEK Optimization Software'});
%Problem Type List
PTYPES = {'sle','lp','milp','bilp','qp','qcqp','miqp','miqcqp','sdp','misdp','snle','scnle','nls','uno','nlp','minlp'};
        
%Default Input Args
if(nargin < 2), err = 1; end %default to create an error
if(nargin < 1), instr = []; end
%Default Internal Args
ret = [];

if(isempty(instr))
    printSolvers();     %print a list of solvers and availability
else
    %Check if we have a known solver
    allS = getAllSolvers(); 
    for i = 1:length(allS)
        if(strcmpi(instr,allS{i}))                    
            ret = check(instr,err);
            if(nargout == 0 && ret), printSolverInfo(instr); end
            return;
        end
    end

    %If we got here, check for a known problem type
    for i = 1:length(PTYPES)
        if(strcmpi(instr,PTYPES{i}))
            ret = findAllForProb(instr,nargout);
            return;
        end
    end
    
    %Still not found, check standard options and alternative names for solvers
    switch(lower(instr))
        %Options
        case 'matrix'
            printSolverMatrix;
        case 'config'
            printConfigMatrix;
        case 'ptypes'
            ret = PTYPES;
        case 'ver'
            ret = genSolverList();
        %Alternative Solver Names
        case 'auto'
            ret = 1; %assume we can find a solver                            
        case 'optim'
            ret = check('matlab',err);
            if(nargout == 0 && ret), printSolverInfo('matlab'); end
        case {'mldivide','fminsearch'}
            ret = 1; %base matlab
        case 'l-bfgs-b'
            ret = check('lbfgsb',err);
            if(nargout == 0 && ret), printSolverInfo('lbfgsb'); end
        case 'lpsolve'
            ret = check('lp_solve',err);
            if(nargout == 0 && ret), printSolverInfo('lp_solve'); end
        case {'lmdif','lm_der','lm_dif'}
            ret = check('lmder',err);
            if(nargout == 0 && ret), printSolverInfo('lmder'); end
        case 'nl2sno'
            ret = check('nl2sol',err);
            if(nargout == 0 && ret), printSolverInfo('nl2sol'); end
        case 'hybrd'
            ret = check('hybrj',err);
            if(nargout == 0 && ret), printSolverInfo('hybrd'); end
        case 'mkltr'
            ret = check('mkltrnls',err);
            if(nargout == 0 && ret), printSolverInfo('mkltrnls'); end
        case 'baron'
            ret = check('baron',err);

        %Check for all, best
        otherwise
            if(strfind(instr,'best'))
                ret = findBestSolver(instr);                        
            elseif(strfind(instr,'all'))
                ret = findAllForProb(instr,nargout);
            else
                error('Unknown Solver or Option: %s',instr); %not found
            end
    end
end

function ok = check(name,err,ptype)
%Check a particular solver is installed
if(nargin < 3), ptype = []; end
switch(name)
    case 'MATLAB'
        if(~isempty(ptype) && strcmpi(ptype,'milp'))                        
            c = which('intlinprog.m'); %only available >= 2014a (8.3)
        elseif(~isempty(ptype) && strcmpi(ptype,'bilp'))
            %bintprog removed in 2014b, used intlinprog instead if available
            if(verLessThan('matlab','8.4')) %2014b
                c = which('bintprog.m');
            else
                c = which('intlinprog.m');
            end                                         
        else
            c = which('fmincon'); %using ver is quite slow
        end
        if(isempty(c))
            if(err)
                error('The MATLAB Optimization Toolbox must be installed');
            else
                ok = 0;
            end
        else
            ok = 1;
        end
        return
    case {'GMATLAB','gmatlab'}
        c = which('patternsearch'); %using ver is quite slow
        if(isempty(c))
            if(err)
                error('The MATLAB Global Optimization Toolbox must be installed');
            else
                ok = 0;
            end
        else
            ok = 1;
        end
        return
    case {'matlab'} %always installed (mldivide)
        ok = 1;
        return
    case {'CPLEX','cplex'} %check cplexlink mex file
        c = which(name);
        %Problems with OPTI + R2013a x64 + Cplex >= v12.5.0.0
        if(~isempty(c) && ~isempty(strfind(computer,'64')))
            if(~isempty(which(['cplexlink1250.' mexext])) || ~isempty(which(['cplexlink1251.' mexext])) || ~isempty(which(['cplexlink1260.' mexext])))
            c = []; %force not available
            end
        end   
    case {'MOSEK','mosek'} %find mex file
        c = which(['mosekopt.' mexext]);
    case {'sedumi','SEDUMI','SeDuMi'}
        %Check matlab file
        c = which('sedumi');
        %Check mex file exists as well
        if(isempty(which(['partitA.' mexext])))
            c = [];
        %Check for v1.31 or higher
        elseif(~isempty(c))
            [~,major,minor] = getVerFromVerFile('sedumi');
            if(major == 1 && minor < 31)
                c = [];
            elseif(major < 1)
                c = [];
            end
        end
    case {'BARON','baron'}
        c = which('baron.m');
        if(isempty(c))
            c = which('baron.p');
        end
    case {'LIPSOL','lipsol'}
        c = which(['ls_blkslv.' mexext]); 
    otherwise
        c = which([name '.' mexext]); %bug fix 22/2/12
        if(~isempty(c)) %ensure mex file is OK
            try
                a = eval(lower(name)); %#ok<NASGU>
            catch
                c = [];
            end
        end
end

if(isempty(c))
    if(err)
        throwAsCaller(MException('opti:checkSolver','OPTI cannot find a valid version of %s. Ensure it is installed on your system.',name));
    else
        ok = 0;
    end
else
    ok = 1;
end

function str = findBestSolver(prob)
%Find the best available solver for a given problem type

global SLE LP MILP BILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP PTYPES %#ok<NUSED>

s = regexp(prob,'_','split');
if(length(s) < 2)
    error('When requesting the best solver it should be of the form best_xxx, e.g. best_lp');
end

%Remove 'D' from dynamic problems
if(s{2}(1) == 'D')
    s{2} = s{2}(2:end);
end

fnd = 0;
for i = 1:length(PTYPES)
    if(strcmpi(s{2},PTYPES{i}))
        try
            list = eval(upper(s{2}));
            fnd = 1;
        catch %#ok<CTCH>
            fnd = 0; %global var doesn't exist?
        end        
        break;
    end
end
%If not found, we don't have a solver for this problem
if(~fnd)
    error('Could not find a problem type match for type %s',upper(s{2}));
%Otherwise return the best solver available
else
    str = [];
    for i = 1:length(list)
        if(check(list{i},0))
            str = list{i};
            break;
        end
    end
    if(isempty(str))
        error('You do not have a solver which can solve a %s explicitly',upper(s{2}));
    end
end

function list = findAllForProb(prob,nout)
%Find all available solvers for a given problem type

global SLE LP MILP BILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP PTYPES %#ok<NUSED>

%Backwards compatibility check
if(strfind(prob,'all_'))
    s = regexp(prob,'_','split');
    prob = s{2};
end

fnd = 0;
for i = 1:length(PTYPES)
    if(strcmpi(prob,PTYPES{i}))
        try
            plist = eval(upper(prob));
            fnd = 1;
        catch %#ok<CTCH>
            fnd = 0; %global var doesn't exist?
        end        
        break;
    end
end
%If found keep availale solvers
if(fnd)
    j = 1; list = [];
    for i = 1:length(plist)
        if(check(plist{i},0,prob))
            list{j} = plist{i}; %#ok<AGROW>
            j = j + 1;
        end
    end
%Check if we want all...    
else
    if(strcmpi(prob,'all'))
        list = getAllSolvers();
        return;
    else
        error('Unknown problem type: %s\nNote the API has changed in v1.34 to just the problem type, e.g. checkSolver(''lp'')',prob);
    end
end

%If not returning anything, just print
if(~nout)
    printSelSolvers(list,upper(prob))
    list = [];
end

function list = getAllSolvers()
%Create a list of all solvers OPTI is interfaced to (not neccesarily found however).
global SLE LP MILP BILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP

fulllist = [SLE LP MILP BILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP];
%Remove duplicates
fulllist = unique(lower(fulllist));
%Sort & Return
list = sort(upper(fulllist));

function printSolvers()
%Print each solver and if it exists

fprintf('\n------------------------------------------------\n')
fprintf('OPTI AVAILABLE SOLVERS:\n\n');
list = genSolverList();
for i = 1:length(list)
    fprintf('%s\n',list{i});
end

fprintf('------------------------------------------------\n')

function list = genSolverList()
%Generate a display compatible list of all solvers, their availability and version number
solv = getAllSolvers();
list = cell(length(solv),1);
for i = 1:length(solv)
    [avail,solVer] = getSolverVer(solv{i});
    list{i} = sprintf('%-10s       %-13s  %s',[solv{i} ':'],avail,solVer);
end

function printSelSolvers(solv,type)
%Print selected solvers

fprintf('\n------------------------------------------------\n')
fprintf(['OPTI ' upper(type) ' SOLVERS:\n']);
if(isempty(solv))
    fprintf('None Available\n');
else
    [~,hd] = getSolverInfo(solv{1},type);
    len = length(strtrim(hd));
    fprintf(' %s\n',hd); j = 1;
    for i = 1:length(solv)
        if(check(solv{i},0)) %only print if the solver is available
            fprintf(['(%2d) %-10s  %-'  num2str(len+4) 's %s\n'],j,upper(solv{i}),getSolverInfo(solv{i},type),getDesc(solv{i},type));
            j = j + 1;
        end
    end
end
fprintf('------------------------------------------------\n')


function printSolverMatrix()
%Print solver vs problem matrix

global SLE LP BILP MILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP PTYPES

fprintf('\n------------------------------------------------\n')
fprintf('OPTI SOLVER vs PROBLEM MATRIX:\n\n');
fprintf(' %-8s  %-5s %-4s %-5s %-6s %-4s %-5s %-5s %-4s %-4s %-4s %-5s %-6s %-5s %-5s %-4s %-4s\n','Solver','SLE','LP','BILP','MILP','QP','QCQP','MIQP','MIQCQP','SDP','MISDP','SNLE','SCNLE','NLS','UNO','NLP','MINLP');

solv = getAllSolvers();
list = {SLE LP BILP MILP QP QCQP MIQP MIQCQP SDP MISDP SNLE SCNLE NLS UNO NLP MINLP};

%REVERSE LOCATION OF BILP AND MILP TO SUIT ABOVE
probs = PTYPES;
idxm = strcmpi(probs,'milp');
idxb = strcmpi(probs,'bilp');
probs{idxm} = 'bilp';
probs{idxb} = 'milp';

for i = 1:length(solv)
    fprintf('%-8s |',solv{i});
    sok = optiSolver(solv{i},0);    
    for j = 1:length(list)
        if(any(strcmpi(solv{i},list{j})))
            if(strcmpi(solv{i},'matlab')) %matlab has optimization toolbox + base packages
                msok = check('MATLAB',0,probs{j});
                if(msok || any(strcmp({'matlab'},list{j})))
                    fprintf('  x  |');
                else
                    fprintf('  o  |');
                end
            else
                if(sok)
                    fprintf('  x  |');
                else
                    fprintf('  o  |');
                end
            end
        else
            fprintf('     |');
        end
    end
    fprintf('\n');
end

fprintf('------------------------------------------------\n')


function printConfigMatrix()
%Print solver vs config matrix

%[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf]
fprintf('\n------------------------------------------------\n')
fprintf('OPTI SOLVER vs CONFIG MATRIX:\n\n');
fprintf('   %-11s %-7s %-7s %-7s %-7s %-7s %-6s %-8s %-7s %-7s %-6s %-7s\n','Solver','miter','meval','mnode','mtime','tolr','tola','tolint','disp','sOpts','ctrlc','iterfnc');

solv = getAllSolvers();

for i = 1:length(solv)
    info = optiSolverInfo(solv{i},'');
    %Separate out variables we use
    miter   = getMarker(info.opt.miter);
    meval   = getMarker(info.opt.meval);
    mnode   = getMarker(info.opt.mnode);
    mtime   = getMarker(info.opt.mtime);
    tolr    = getMarker(info.opt.tolr);
    tola    = getMarker(info.opt.tola);
    tolint  = getMarker(info.opt.tolint);
    disp    = getMarker(info.opt.disp);
    sOpts   = getMarker(info.opt.opts);
    ctrlc   = getMarker(info.opt.ctrlc);
    iterf   = getMarker(info.opt.iterf);    
    
    %Print line
    fprintf('%-11s  |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |   %1s   |\n',solv{i},miter,meval,mnode,mtime,tolr,tola,tolint,disp,sOpts,ctrlc,iterf);
end

fprintf('------------------------------------------------\n')


function str = getDesc(solv,type)
%Return solver description
global DESC

str = DESC(lower(solv));
if(strcmp(solv,'MATLAB'))   
    str = [str getMatlabDesc(type)];
elseif(strcmp(solv,'GMATLAB'))   
    str = [str getGMatlabDesc(type)];    
elseif(strcmp(solv,'matlab'))
    str = getMatlabDesc(type);
end

function str = getMatlabDesc(type)
%Return MATLAB solver

switch(lower(type))
    case 'sle'
        str = 'mldivide';
    case 'lp'
        str = 'linprog';
    case 'qp'
        str = 'quadprog';
    case 'bilp'
        if(verLessThan('matlab','8.4')) %2014b
            str = 'bintprog';
        else
            str = 'intlinprog';
        end
    case 'milp'
        str = 'intlinprog';
    case 'nlp'
        str = 'fmincon';
    case 'snle'
        str = 'fsolve';
    case 'nls'
        str = 'lsqnonlin';
    case 'uno'
        str = 'fminsearch';
    otherwise
        str = '';
end

function str = getGMatlabDesc(type)
%Return Global MATLAB Solver

switch(lower(type))
    case {'uno','nlp'}
        str = 'patternsearch';
    case 'minlp'
        str = 'ga';
    otherwise
        str = '';
end
        
function [avail,solVer] = getSolverVer(solver)
%Return currently installed solver version

ok = check(solver,0);
if(ok)
    avail = 'Available';    
    switch(lower(solver))
        %Interfaced solvers
        case 'cplex'
            c = Cplex;
            solVer = ['v' c.getVersion];
            clear c;
        case 'matlab'
            c = ver('optim');
            if(isfield(c,'Version'))
                solVer = ['v' c.Version];
            else
                solVer = '';
            end
            clear c;
        case 'gmatlab'
            c = ver('globaloptim');
            if(isfield(c,'Version'))
                solVer = ['v' c.Version];
            else
                solVer = '';
            end
            clear c;
        case 'mosek'
            solVer = ''; %unknown for the moment
        case 'sedumi'
            solVer = getVerFromVerFile('sedumi');
        case 'lipsol'
            solVer = getVerFromVerFile('lipsol');
        case 'opti'
            c = optiver;
            solVer = ['v' num2str(c)];  
        case 'baron' %unknown
            solVer = '';
        %Supplied Solvers               
        otherwise
            try
                v = feval(lower(solver));
                if(~strcmp(v,''))
                    solVer = ['v' v];
                else
                    solVer = '';
                end
            catch %#ok<CTCH>
                solVer = '';
            end
    end        
else
    avail = 'Not Available';    
    solVer = '';
end


function [str,hd] = getSolverInfo(solver,type)
%Generate a string of solver information for a particular solver

%Get Solver Information
info = optiSolverInfo(solver,type);
%Separate out variables we use
bnd     = getMarker(info.con.bnd);
lineq   = getMarker(info.con.lineq);
leq     = getMarker(info.con.leq);
qineq   = getMarker(info.con.qineq);
ss   = getMarker(info.con.sos);
sdp     = getMarker(info.con.sdp);
nlineq  = getMarker(info.con.nlineq);
nleq    = getMarker(info.con.nleq);
der1    = getMarker(info.der1);
der2    = getMarker(info.der2);
glb     = getMarker(info.global);
sp      = getMarker(info.sparse);

lhs = '%-17s';
switch(lower(type))
    case {'lp','qp'}
        hd = sprintf([lhs 'BD LI LE SP'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s',bnd,lineq,leq,sp);
    case {'milp','miqp'}
        hd = sprintf([lhs 'BD LI LE SS SP'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s %-2s',bnd,lineq,leq,ss,sp);
    case 'bilp'
        hd = sprintf([lhs 'LI LE SP'],'');
        str = sprintf(' %-2s %-2s %-2s',lineq,leq,sp);        
    case {'qcqp','miqcqp'}
        hd = sprintf([lhs 'BD LI LE QI SP'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s',bnd,lineq,leq,qineq,sp);
    case {'sdp','misdp'}
        hd = sprintf([lhs 'BD LI LE QI SD SP'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s %-2s',bnd,lineq,leq,qineq,sdp,sp);
    case 'nls'
        hd = sprintf([lhs 'BD LI LE D1'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s',bnd,lineq,leq,der1);
    case {'nlp','minlp'}
        hd = sprintf([lhs 'BD LI LE QI NI NE SP D1 D2 GL'],'');
        str = sprintf(' %-2s %-2s %-2s %-2s %-2s %-2s %-2s',bnd,lineq,leq,qineq,nlineq,nleq,sp,der1,der2,glb);
    otherwise %no constraints
        hd = ''; 
        str = '';
end    

function str = getMarker(value)
switch(value)
    case 1
        str = 'x';
    case 0
        str = '';
    case -1
        str = '*';
end

%Print Solver Information
function printSolverInfo(solver)
fprintf('\n-----------------------------------------------------------\n');
%Print heading
[~,solVer] = getSolverVer(solver);
if(strcmpi(solver,'MATLAB'))
    solverstr = 'MATLAB [Optimization Toolbox]';
elseif(strcmpi(solver,'GMATLAB'))
    solverstr = 'GMATLAB [Global Optimization Toolbox]';
else
    solverstr = upper(solver);
end
fprintf('OPTI SOLVER: %s [%s]\n',solverstr,solVer);  
fprintf('-----------------------------------------------------------\n');
%Get solver info
info = optiSolverInfo(solver);
%Print Special Info
if(info.der1 || info.der2 || info.sparse || info.global || info.multialg)
   fprintf('Solver Features:\n');
   if(info.parallel == 1), fprintf('\t- Multi-Threaded (Parallel) Solver\n'); end
   if(info.multialg == 1), fprintf('\t- Contains Multiple Algorithms\n'); end 
   if(info.whitebox == 1), fprintf('\t- White-Box NL Solver Which Exploits Problem Structure\n'); end 
   if(info.global == 1)
       fprintf('\t- Solves Global Optimization Problems\n');
   elseif(info.global == -1)
       fprintf('\t- Some Algorithms Solve Global Optimization Problems\n'); 
   end          
   if(info.con.int), fprintf('\t- Solves Mixed Integer Problems\n'); end
   if(info.sparse), fprintf('\t- Solves Sparse Problems\n'); end
   if(info.der1 == 1)
       fprintf('\t- Requires First Derivatives\n'); 
   elseif(info.der1 == -1)
       if(info.multialg==1)
            fprintf('\t- Some Algorithms Require First Derivatives\n'); 
       else
            fprintf('\t- Can Accept First Derivatives (Or Solver Will Approximate)\n');
       end
   end
   if(info.der2 == 1)
       fprintf('\t- Requires Second Derivatives\n'); 
   elseif(info.der2 == -1)
       fprintf('\t- Can Accept Second Derivatives (Or Solver Will Approximate)\n');
   end
   fprintf('-----------------------------------------------------------\n'); 
end
%Print Constraint Info
if(info.con.con)
    fprintf('Supported Constraints:\n');
    printCon('Decision Variable Bounds',info.con.bnd);
    printCon('Linear Inequalities',info.con.lineq);
    printCon('Linear Equalities',info.con.leq);
    printCon('Quadratic Inequalities',info.con.qineq);
    printCon('Quadratic Equalities',info.con.qeq);
    printCon('Semidefinite (Linear Matrix Inequalities)',info.con.sdp);
    printCon('Nonlinear Inequalities',info.con.nlineq);
    printCon('Nonlinear Equalities',info.con.nleq);
    printCon('Integer Variables',info.con.int);
    printCon('Special Ordered Sets',info.con.sos);
    fprintf('-----------------------------------------------------------\n');
end
%See if we have any unique info
u=struct2cell(info.unq); uu = 0;
for i = 1:length(u)
    if(ischar(u{i}) || any(u{i} ~= 0))
        uu = 1;
        break;
    end
end
%Print Unique Info
if(uu)
    fprintf('Special Requirements / Options:\n');
    if(info.unq.x0 == 1)
        fprintf('\t- Requires an initial solution guess via x0\n');
    elseif(info.unq.x0 == -1)
        fprintf('\t- Can accept an initial solution guess via x0\n');
    end
    switch(info.unq.qpform)
        case 'sym'
            fprintf('\t- QP H Matrix must be symmetric\n');
        case 'tril'
            fprintf('\t- QP H Matrix must be symmetric lower triangular\n');
        case 'triu'
            fprintf('\t- QP H Matrix must be symmetric upper triangular\n');
    end
    switch(info.unq.hform)
        case 'sym'
            fprintf('\t- Nonlinear Hessian (if supplied) must be symmetric\n');
        case 'tril'
            fprintf('\t- Nonlinear Hessian (if supplied) must be symmetric lower triangular\n');
        case 'triu'
            fprintf('\t- Nonlinear Hessian (if supplied) must be symmetric upper triangular\n');
    end
    if(info.unq.jstr), fprintf('\t- Requires the structure of the Jacobian\n'); end       
    if(info.unq.hstr), fprintf('\t- Requires the structure of the Hessian (if supplied)\n'); end       
    fprintf('-----------------------------------------------------------\n');
end
        
%Print info from MEX file
switch(lower(solver))
    case {'cplex','mosek','matlab','gmatlab','sedumi','lipsol','baron'}
        %Nothing
    case {'filtersd'}
        fprintf('MEX FILE DETAILS:');
        eval([lower(solver) 'sp']);
    otherwise
        fprintf('MEX FILE DETAILS:');
        eval(lower(solver));
end

function printCon(con,val)
if(val == 1)
    fprintf('\t- %s\n',con); 
elseif(val == -1)
    fprintf('\t- %s (but may be converted internally)\n',con);
end

function [solVer,major,minor] = getVerFromVerFile(solverP)
solVer = ''; major = 0; minor = 0;
%Find sedumi.m and Version .txt
p = which(solverP); p2 = which('Version.txt','-all');
if(~isempty(p) && ~isempty(p2))
    %Find final file sep in solver path
    ind = strfind(p,filesep);
    %For each version path, see if we can find a matching solver path
    for i = 1:length(p2)
        ind2 = strfind(p2{i},filesep);
        %compare paths up to find file sep
        if(strcmp(p(1:ind(end)),p2{i}(1:ind2(end))))
            %open and reading version
            fid = fopen(p2{i});
            if(fid > 0)
                v = fscanf(fid,'%d.%d');
                fclose(fid);
                major = v(1); minor = v(2);
                solVer = sprintf('v%d.%d',v(1),v(2));
            else
                fclose(fid);
            end
            break;
        end
    end
end

