function info = optiSolverInfo(solver,probtype,nlprob,opts)
%OPTISOLVERINFO  Return Configuration Information about a particular solver
%
%   info = optiSolverInfo(solver,probtype)

%Default input args and checks
if(nargin < 4), opts = []; end
if(nargin < 3), nlprob = []; end
if(nargin < 2 || isempty(probtype)), probtype = ''; end
if(nargin < 1 || ~ischar(solver)), error('You must supply the name of the solver to this function'); end

%Create default information structure
info.der1 = 0;
info.der2 = 0;
info.global = 0;
info.sparse = 0;
info.multialg = 0;
info.parallel = 0;
info.whitebox = 0;
info.con = struct('con',0,'bnd',0,'lineq',0,'leq',0,'qineq',0,'qeq',0,'sdp',0,'nlineq',0,...
                  'nleq',0,'int',0,'sos',0);
info.opt = struct('miter',0,'meval',0,'mnode',0,'mtime',0,'tolr',0,...
                  'tola',0,'tolint',0,'disp',0,'opts',0,'ctrlc',0,...
                  'iterf',0);       
info.unq = struct('x0',0,'qpform',0,'hform',0,'jstr',0,'hstr',0);
              
%Constants
YES = 1;
SOME = -1;
              
%Fill in information based on selected solver
switch(lower(solver))
    case 'baron'
        info.global = YES;
        info.sparse = YES;
        info.whitebox = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,yes,no,yes,yes,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
%         info.con.qineq = YES;
        info.con.nlineq = YES;
        info.con.nleq = YES;
        info.con.int = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(some,no,yes,yes,yes,no,no,yes,no,yes,no);
        info.opt.miter = SOME;
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tola = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = SOME;
    
    case 'bonmin'
        info.der1 = YES;
        info.der2 = SOME;
        info.sparse = YES;
        info.multialg = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,some,some,no,no,yes,yes,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.nlineq = YES;
        info.con.nleq = YES;     
        info.con.int = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(some,no,yes,yes,some,no,yes,yes,yes,yes,no);
        info.opt.miter = SOME;
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tolr = SOME;
        info.opt.tolint = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.x0 = YES;
        info.unq.hform = 'tril';
        info.unq.jstr = YES;
        info.unq.hstr = YES;
        
    case 'clp'
        info.sparse = YES;
        info.multialg = YES;
        info.parallel = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);  
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,no,yes,no,yes,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp =  YES;
        info.opt.ctrlc = YES;
        info.opt.opts = YES;
        %unique
        info.unq.qpform = 'tril';
        
    case 'cbc'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);   
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.int = YES;
        info.con.sos = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,no,yes,yes,no,no,yes,yes,no,yes,no);
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tolint = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.opts = YES;
        %unique
        info.unq.x0 = SOME;
        info.unq.qpform = 'tril';
        
    case 'cplex'
        info.sparse = YES;
        info.parallel = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,yes,no,no,no,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.qineq = YES;
        info.con.int = YES;
        info.con.sos = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,yes,yes,yes,no,yes,yes,yes,yes,no);
        info.opt.miter = YES;
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tolint = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.x0 = SOME;
        info.unq.qpform = 'sym';
        
    case 'csdp'
        info.sparse = YES;
        info.parallel = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,some,no,yes,no,no,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = SOME;
        info.con.sdp = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,no,yes,yes,yes,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.x0 = SOME;
        
    case 'dsdp'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,some,no,yes,no,no,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = SOME;
        info.con.sdp = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,no,yes,no,yes,yes,yes,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tola = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.x0 = SOME;
        
    case 'filtersd'
        info.der1 = YES;
        info.sparse = YES;
        info.multialg = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,some,some,no,no,yes,yes,yes); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = SOME;
        info.con.leq = SOME;
        info.con.nlineq = YES;
        info.con.nleq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,no,yes,no,some,no,yes,no,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tola = SOME;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        info.unq.jstr = YES;
        
    case 'glpk'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);  
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.int = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,yes,yes,no,yes,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tolint = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        
    case 'gmatlab'   
        info.global = YES;
        info.sparse = YES;
        info.multialg = YES;
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.nlineq = YES;     
        info.con.int = YES;
        switch(lower(probtype))
            case {'uno','nlp',''}
                %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,yes,yes,yes);                 
                info.con.leq = YES;               
                info.con.nleq = YES;
        end
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,some,yes,yes,no,some,yes,yes,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mnode = SOME;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tolint = SOME;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = SOME;
        
    case 'hybrj'
        info.der1 = SOME;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,yes,no,yes,no,no,no,yes,no,yes,yes);
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;  
        %unique
        info.unq.x0 = YES;
        
    case 'ipopt'
        info.der1 = YES;
        info.der2 = SOME;    
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,some,some,no,no,yes,yes,yes);   
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.nlineq = YES;
        info.con.nleq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,no,yes,yes,yes,yes);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        info.unq.hform = 'tril';
        info.unq.jstr = YES;
        info.unq.hstr = YES;
        
    case 'lbfgsb'
        info.der1 = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,no,no,no,no,no,no,no); 
        info.con.con = YES;
        info.con.bnd = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,no,yes,no,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        info.opt.opts = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'levmar'
        info.der1 = YES;        
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,no);    
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,no,no,no,no,yes,no,some,yes);
        info.opt.miter = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = SOME;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'lipsol'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);  
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = SOME;
        info.con.leq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,no,no,yes,no,yes,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp =  YES;
        info.opt.ctrlc = YES;
        
    case 'lmder'
        info.der1 = SOME;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,yes,no,yes,yes,no,no,yes,no,yes,yes);
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'lp_solve'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;     
        info.con.int = YES;
        info.con.sos = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,no,some,yes,no,no,yes,yes,no,no,no);
        info.opt.mnode = SOME;
        info.opt.mtime = YES;
        info.opt.tolint = YES;
        info.opt.disp = YES;       
        
    case 'm1qn3'
        info.der1 = YES;        
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,no,yes,no,some,no,yes,no,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tola = SOME;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'matlab'
        info.der1 = SOME;
        info.der2 = SOME;
        info.sparse = YES;
        info.multialg = YES;
        info.con.con = YES;
        info.con.bnd = YES;                 
        switch(lower(probtype))
            case {'lp','qp','bilp','milp'}
                %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);
                info.con.lineq = YES;
                info.con.leq = YES;
                %unique
                info.unq.x0 = SOME;
                if(~isempty(which('intlinprog.m')))
                    info.con.int = YES;
                end
            case {'nlp',''}
                info.der1 = YES;
                info.der2 = SOME;
                %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,yes,yes,yes);                 
                info.con.lineq = YES;
                info.con.leq = YES;
                info.con.nlineq = YES;
                info.con.nleq = YES;
                if(~isempty(which('intlinprog.m')))
                    info.con.int = YES;
                end
                %unique
                info.unq.x0 = YES;
                info.unq.hform = 'sym';
                info.unq.jstr = YES;
                info.unq.hstr = YES;
        end
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,some,yes,yes,no,some,yes,yes,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        
    case 'mkltrnls'
        %If we don't have options, assume not called by OPTI and don't allow internal mklJac
        if(isempty(opts))
           info.der1 = SOME;
        end
        %bnd = yes;
        info.con.con = YES;
        info.con.bnd = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,yes,yes,yes,no,yes,no,yes,yes); 
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tola = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'mosek'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,yes,no,no,no,yes); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.qineq = YES;    
        info.con.int = YES;
        info.con.sos = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,yes,yes,yes,no,yes,yes,yes,yes,no);
        info.opt.miter = YES;
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tolint = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.x0 = SOME;
        info.unq.qpform = 'tril';
        
    case 'mumps'
        %Nothing!
        
    case 'nl2sol'
        info.der1 = SOME;
        %bnd = yes;
        info.con.con = YES;
        info.con.bnd = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,no,no,yes,yes,no,yes,no,some,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.tolr = YES;
        info.opt.tola = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = SOME;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'nlopt'
        %Some NLOPT algorithms don't require derivatives, see if we have problem data
        if(~isempty(opts))       
            if(isfield(opts,'solverOpts') && ~isempty(opts.solverOpts) && isfield(opts.solverOpts,'algorithm'))
                [~,~,~,~,~,info.der1] = nloptSolver(opts.solverOpts.algorithm);
            elseif(~isempty(nlprob) && isfield(nlprob,'algorithm'))
                [~,~,~,~,~,info.der1] = nloptSolver(nlprob.algorithm);
            end %else assume not required (default based on solver)
        else
            info.der1 = SOME; %generally some do
        end
        info.global = SOME;
        info.multialg = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,some,some,no,no,yes,yes,no);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = SOME;
        info.con.leq = SOME;
        info.con.nlineq = YES;
        info.con.nleq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,yes,no,yes,yes,yes,no,yes,yes,some,yes);
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.tola = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = SOME;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'nomad'
        info.global = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,some,no,no,no,yes,no,no); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = SOME;
        info.con.nlineq = YES;
        info.con.int = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,no,yes,yes,no,no,yes,yes,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'ooqp'
        info.sparse = YES;
        info.multialg = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(no,no,no,no,no,no,no,yes,no,no,no);
        info.opt.miter = YES;
        info.opt.mtime = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        %unique
        info.unq.qpform = 'tril';
        
    case 'pswarm'
        info.global = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,no,no,no,no,no,no); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,yes,no,yes,yes,no,no,yes,yes,yes,yes);
        info.opt.miter = YES;
        info.opt.meval = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.disp = YES;
        info.opt.opts = YES;
        info.opt.ctrlc = YES;
        info.opt.iterf = YES;
        %unique
        info.unq.x0 = YES;
        
    case 'qsopt'
        info.sparse = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,no,no,no,no,yes); 
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(yes,no,no,some,no,no,no,some,no,no,no);
        info.opt.miter = YES;
        info.opt.mtime = SOME;
        info.opt.disp = SOME;        
        
    case 'scip'
        info.global = YES;
        info.sparse = YES;
        info.whitebox = YES;
        %[bnd,lineq,leq,qineq,sdp,nlineq,nleq,sp] = fillFields(yes,yes,yes,yes,no,yes,yes,yes);
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = YES;
        info.con.qineq = YES;
        info.con.qeq = YES;
        info.con.nlineq = YES;
        info.con.nleq = YES;
        info.con.int = YES;
        info.con.sos = YES;
        %[miter,meval,mnode,mtime,tolr,tola,tolint,disp,opts,ctrlc,iterf] = fillCFields(some,no,yes,yes,yes,no,no,yes,no,yes,no);
        info.opt.miter = SOME;
        info.opt.mnode = YES;
        info.opt.mtime = YES;
        info.opt.tolr = YES;
        info.opt.opts = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        
    case 'sedumi'
        info.sparse = YES;
        info.multialg = YES;
        %Constraints
        info.con.con = YES;
        info.con.bnd = YES;
        info.con.lineq = YES;
        info.con.leq = SOME;
        info.con.sdp = YES;
        %Options
        info.opt.miter = YES;
        info.opt.disp = YES;
        info.opt.ctrlc = YES;
        
    case 'auto'
        %Assume we require derivatives
        info.der1 = YES;
        
    otherwise
        error('Unknown solver: %s',solver);
end
              