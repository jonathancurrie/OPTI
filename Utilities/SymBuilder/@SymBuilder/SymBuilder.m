classdef SymBuilder < handle
%SYMBUILDER  Create a SYMBUILDER object for creating symbolic optimization problems
%
%   B = SymBuilder() creates a blank optimization problem
%
%   See also SymBuilder.AddObj SymBuilder.AddCon SymBuilder.Build
%
%   Copyright (C) 2012 Jonathan Currie (www.inverseproblem.co.nz)
    
    properties(SetAccess=private)
        vars        %Symbolic array of symbolic variables
        sobj        %Symbolic matrix of equations
        jac         %Symbolic jacobian of equations
        hess        %Structure of hessians of each equation
        hesslag     %Symbolic hessian of the lagrangian  
        Opt         %OPTI object
    end
    
    properties(SetAccess=private)%,GetAccess=private)
        eqs         %Equations symbolic array
        cl          %Constraint Lower Bounds
        cu          %Constraint Upper Bounds
        constnt     %Cell array of constants
        exprsn      %Cell array of expressions
        bnds        %Cell array of bounds
        lb          %Double array of lower bounds
        ub          %Double array of upper bounds
        vartypes    %Cell array of variable types
        xtype       %Char array of variable types
        indobj      %Index of objectives in equations
        objlin      %Index of linear, quadratic and nonlinear objectives
        conlin      %Index of linear, quadratic and nonlinear constraints
        objbias     %Problem bias terms (for linear and quadratic problems)
        noObjs      %Counter of number objectives
        noCons      %Counter of number of constraints
        bldstat     %Build Status
        bldtime     %Build Time
        resgrp      %Cell array of result groups
        resexp      %Cell array of result expressions
        verbose     %Verbosity mode
    end
    
    methods
        function B = SymBuilder(verbose)
            if(nargin < 1), verbose = true; end
            if (~exist('syms.m','file')), error('The Symbolic Math Toolbox must be installed to use SymBuilder!'); end
            %Initialization
            B.noObjs = 0;
            B.noCons = 0;
            B.verbose = verbose;
            Unbuild(B);            
        end
                
        function B = AddObj(B,str)
        %Add a General Objective (Linear or Nonlinear)
        
            %Check for string
            if(ischar(str))
                digits(16);
                wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
                s = sym(str);
                warning(wstate);
            elseif(isa(str,'sym'))
                s = str;
            else
                error('Unknown objective form');
            end
        
            %Concatenate with existing equations
            B.eqs = [B.eqs;s];
            B.cl = [B.cl;NaN];
            B.cu = [B.cu;NaN];
            if(isempty(B.indobj))
                B.indobj = length(B.cl);
            else
                B.indobj = [B.indobj;length(B.cl)];
            end
            %Add no objs
            B.noObjs = B.noObjs + 1;  
            %Indicate Rebuild Required
            Unbuild(B);
        end
                
        function B = AddCon(B,str_in,cl,cu)
        %Add a General Constraint (Linear or Nonlinear)
        
            %Check for string
            if(ischar(str_in))
                %Parse constraint string
                [str,l,u] = SymBuilder.parseConstraint(str_in);
                digits(16);
                wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
                s = sym(str);
                warning(wstate);
            elseif(isa(str_in,'sym'))
                if(nargin < 4 || ~isnumeric(cl) || ~isnumeric(cu) || length(cl)~=length(cu) || length(cl)~=length(str_in))
                    error('When supplying a symbolic expression to AddCon, you must supply as AddCon(sym_con,cl,cu) with cl, cu as numeric vectors of the same length.');
                end
                if(size(str_in,2) > 1)
                    s = str_in.';
                else
                    s = str_in;
                end
                l = cl; u = cu;               
            else        
                 error('Unknown constraint form');    
            end
            %Concatenate with existing equations if not empty
            if(~isempty(s))
                B.eqs = [B.eqs;s];
                B.cl = [B.cl;l];
                B.cu = [B.cu;u];           
                %Add no cons
                B.noCons = B.noCons + length(l);  
                %Indicate Rebuild Required
                Unbuild(B);
            else
                optiwarn('symb:emptycon','The following constraint will be ignored by SymBuilder as both sides appear to be numeric.\n  %s\n',str_in);
            end
        end
        
        function B = AddLinCon(B,x,A,rl,ru)
        %Add a Linear Constraint
        
            if(nargin < 5 || ~isnumeric(rl) || ~isnumeric(ru) || length(rl)~=length(ru) || length(rl)~=size(A,1))
                error('When supplying linear constraints, you must supply as AddLinCon(x,A,rl,ru) with rl, ru as numeric vectors of the same length.');
            elseif(~isa(x,'sym') || length(x) ~= size(A,2))
                error('x must be a symbolic variable vector = size(A,2)');
            end        
            %Evaluate Matrix & Store            
            B.eqs = [B.eqs;A*x];
            B.cl = [B.cl;rl];
            B.cu = [B.cu;ru];
            %Add no cons
            B.noCons = B.noCons + length(rl);  
            %Indicate Rebuild Required
            Unbuild(B);
        end
        
        
        function B = AddConstraint(B,str)
        %Same as AddCon()
            B = AddCon(B,str);
        end
                
        function B = AddConstant(B,str,val,varargin)
        %Identify variable as a constant 
        
            %Concatenate with existing constants
            B.constnt = [B.constnt; {str val}];   
            %If user has passed multiple constants, process all
            if(nargin > 3 && ~isempty(varargin))
                if(mod(length(varargin),2))
                    error('You must supply a name and value for each constant');
                else
                    for i = 1:2:length(varargin)
                        B.constnt = [B.constnt; varargin(i:i+1)];
                    end
                end
            end
            %Indicate Rebuild Required
            Unbuild(B);
        end  
                
        function B = AddExpression(B,str,grp,ex)
        %Identify variable as an expression
        
            %Parse expression string
            [str,exp] = SymBuilder.parseExpression(str);
            %Concatenate with existing expressions
            B.exprsn = [B.exprsn; {str exp}]; 
            %If group passed, add to result group
            if(nargin > 3), B = AddResultExp(B,grp,exp,ex);
            elseif(nargin > 2), B = AddResultExp(B,grp,exp); end
            %Indicate Rebuild Required
            Unbuild(B);
        end
        
        function B = AddBound(B,str)
        %Add variable bound
        
            %Parse expression string
            [var,llb,lub] = SymBuilder.parseBound(str);
            %Concatenate with existing bounds
            B.bnds = [B.bnds; {var llb lub}];             
            %Indicate Rebuild Required
            Unbuild(B);        
        end
        
        function B = AddBounds(B,x,lb,ub)
        %Add numerical bounds
            
            bds = cell(length(x),3);
            for i = 1:length(x)
                bds{i,1} = char(x(i));
                bds{i,2} = sprintf('%1.15g',lb(i));
                bds{i,3} = sprintf('%1.15g',ub(i));
            end        
            %Concatenate with existing bounds
            B.bnds = [B.bnds; bds];             
            %Indicate Rebuild Required
            Unbuild(B); 
        end
            
        function B = AddInteger(B,str)
        %Add integer constraint
        
            %Parse expression string
            [var,xx] = SymBuilder.parseInteger(str);
            %Concatenate with existing constants
            B.vartypes = [B.vartypes; {var xx}];             
            %Indicate Rebuild Required
            Unbuild(B);        
        end

        function B = AddIntegers(B,x,str)
        %Add integer constraints with character string
        
            xx = cell(length(x),2);
            for i = 1:length(x)
                xx{i,1} = char(x(i));
                xx{i,2} = str(i);
            end 
            %Concatenate with existing constants
            B.vartypes = [B.vartypes; xx];             
            %Indicate Rebuild Required
            Unbuild(B);        
        end
        
        function B = AddResultGroup(B,name,str)
        %Add Result Group
            %Concatenate with existing groups
            B.resgrp = [B.resgrp; {name str}];        
        end
        
        function B = AddResultExp(B,name,str,bin)
        %Add Result Expression            
            %Parse expression string
            [group,name] = SymBuilder.parseResExp(name);
            %Concatenate with existing groups
            if(nargin > 3)
                B.resexp = [B.resexp; {group name str bin}];
            else
                B.resexp = [B.resexp; {group name str []}];        
            end
        end
                
        function B = Draft(B)
        %Draft the object by solving for the system Jacobian
            t = tic;
            %Build Symbolic representation of the equations
            if(B.verbose), fprintf('\nGenerating Symbolic Representation of Equations...'); end
            buildSymRep(B);
            if(B.verbose), fprintf('Done\n'); end
            %Generate Equations Jacobian
            if(B.verbose), fprintf('Generating Symbolic Jacobian...'); end
            buildJac(B);
            if(B.verbose), fprintf('Done\n'); end
            %Determine Linear & Nonlinear Equations
            if(B.verbose), fprintf('Generating Equation Linearity...'); end
            detLinearity(B);
            if(B.verbose), fprintf('Done\n'); end
            %Save build time
            B.bldtime = toc(t);
            %Assign build Status
            B.bldstat = 'draft';
        end
                
        function B = Build(B)
        %Build the object by solving for the system Hessian
            t = tic;
            %Build Symbolic representation of the equations
            if(B.verbose), fprintf('\nGenerating Symbolic Representation of Equations...'); end
            buildSymRep(B);
            if(B.verbose), fprintf('Done\n'); end
            %Generate Equations Jacobian
            if(B.verbose), fprintf('Generating Symbolic Jacobian...'); end
            buildJac(B);
            if(B.verbose), fprintf('Done\n'); end
            %Generate Equations Hessians (also determines linearity)
            if(B.verbose), fprintf('Generating Symbolic Hessian...'); end
            buildHess(B);
            if(B.verbose), fprintf('Done\n'); end
            %Save build time
            B.bldtime = toc(t);
            %Assign build Status
            B.bldstat = 'built';
        end
        
        function B = Unbuild(B)
        %Unbuild object (free memory)
        
            B.jac = [];
            B.hess = [];
            B.hesslag = [];
            B.bldstat = 'unbuilt';
            B.bldtime = 0;
            B.Opt = [];
            
            clear symb_cb
        end
              
        function prob = GetLinProb(B)
        %Get Linear Problem Matrices [min f'*x s.t. rl <= A*x <= ru]    
            
            %Check is built
            IsBuilt(B,1); 
            prob = optiprob; 
            %Get Linear Objective Index
            ind = B.objlin == 1;
            if(sum(ind) == 0)
                error('This object does not contain a linear objective');
            end
            prob.f = double(B.jac(ind,:))';
            %Get Linear Constraint Indices
            ind = B.conlin == 1; 
            %Check we have some constraints
            if(B.noCons && sum(ind) > 0)
                if(isempty(B.objbias)), B.objbias = zeros(size(ind)); end %HACK
                %Get Linear Jacobian
                prob.A = sparse(double(B.jac(ind,:)));
                prob.rl = B.cl(ind) - B.objbias(ind);
                prob.ru = B.cu(ind) - B.objbias(ind);
            end
            %Bounds & Integer Vars
            prob.lb = B.lb;
            prob.ub = B.ub;
            prob.int = B.xtype; 
            prob.objbias = B.objbias(B.indobj);
            %Problem Name
            if(any(B.xtype=='I'))
                prob.Name = 'SymBuilder MILP';
            elseif(any(B.xtype=='B'))
                prob.Name = 'SymBuilder BILP';
            else
                prob.Name = 'SymBuilder LP';
            end                       
        end
        
        function prob = GetQuadProb(B)
        %Get Quadratic Problem Matrices [min 0.5*x'*H'x + f'*x 
        %s.t. rl <= A*x <= ru, qrl <= x'*Q*x + l'*x <= qru]
            
            %Check is built
            IsBuilt(B,1); 
            prob = optiprob;
            %Get quadratic objective indices
            ind = B.objlin == 2;
            if(isempty(ind))
                error('This object does not contain an objective');
            elseif(sum(ind) > 1)
                error('This interface only supports a single quadratic objective');
            end
            %Get Hessian (DON'T NEED 0.5!)
            if(~any(ind)) %linear objective
                ind = B.objlin == 1;
                n = size(B.jac,2);
                prob.H = spalloc(n,n,0);
                prob.f = double(B.jac(ind,:)); 
            else
                prob.H = sparse(double(B.hess(ind).H));
                %Get Gradient, remove 2nd derivative parts
                f = B.jac(ind,:); v = symvar(f);
                prob.f = double(subs(f,v,{zeros(size(v))}).');
            end
            
            %Get Linear Constraint Indices
            ind = B.conlin == 1; 
            %Check we have some constraints
            if(B.noCons && sum(ind) > 0)
                %Get Linear Jacobian
                prob.A = sparse(double(B.jac(ind,:)));
                prob.rl = B.cl(ind) - B.objbias(ind);
                prob.ru = B.cu(ind) - B.objbias(ind);
            end
            %Get Quadratic Constraint Indices
            ind = B.conlin == 2; 
            %Check we have some constraints
            if(B.noCons && sum(ind) > 0)
                prob.qrl = B.cl(ind) - B.objbias(ind);
                prob.qru = B.cu(ind) - B.objbias(ind);
                ind = find(ind);
                %Get each quadratic constraint
                len = length(ind);
                prob.Q = cell(len,1);
                prob.l = cell(len,1);
                for i = 1:len
                    %Get Hessian (0.5* requried!)
                    prob.Q{i} = 0.5*sparse(double(B.hess(ind(i)).H));
                    %Get Gradient, remove 2nd derivative parts
                    f = B.jac(ind(i),:); v = symvar(f);
                    prob.l{i} = double(subs(f,v,{zeros(size(v))}).');
                end                
            end            
            %Bounds & Integer Vars
            prob.lb = B.lb;
            prob.ub = B.ub;
            prob.int = B.xtype; 
            prob.objbias = B.objbias(B.indobj);
            %Problem Name
            if(any(B.xtype=='I'))
                prob.Name = 'SymBuilder MIQP';
            else
                prob.Name = 'SymBuilder QP';
            end 
        end
        
        function prob = GetNLProb(B,opts)
        %Get Nonlinear Problem Callbacks, generating files as we go
            
            if(nargin < 2), opts = symbset; end
            opts.verbose = B.verbose;
            %Check is built
            IsBuilt(B,1);
            prob = optiprob;   
            %Determine most efficient callback form IF auto selected
            if(strcmpi(opts.cbmode,'auto'))        
                %Only relevant if we have any nonlinear terms & > 8 variables
                if((any(B.objlin > 2) || any(B.conlin > 2)) && length(B.vars) > 8)
                    %See if we have 2nd deriv info
                    if(~isempty(B.hess))
                        %Generate Hess Lag (will need it later anyway)
                        buildHessLag(B);                  
                        %Build Hessian Structure
                        lHstr = zeros(size(B.hesslag));
                        lHstr(logical(B.hesslag ~= 0)) = 1;
                        %If greater than 10% dense, use CppAD
                        if(nnz(lHstr)/numel(lHstr) > 0.1)
                            if(SymBuilder.CheckCppADCompile(false))
                                opts.cbmode = 'cppad';
                            else
                                opts.cbmode = 'mcode';
                            end
                        else
                            opts.cbmode = 'mcode';
                        end
                    %Only 1st deriv info
                    else
                        ind = B.conlin > 0; 
                        ajac = B.jac(ind,:);
                        %Build Jacobian Structure
                        jacstr = zeros(size(ajac));
                        jacstr(logical(ajac ~= 0)) = 1;
                        %If greater than 10% dense, use CppAD
                        if(nnz(jacstr)/numel(jacstr) > 0.1)
                            if(SymBuilder.CheckCppADCompile(false))
                                opts.cbmode = 'cppad';
                            else
                                opts.cbmode = 'mcode';
                            end
                        else
                            opts.cbmode = 'mcode';
                        end
                    end
                else
                    opts.cbmode = 'mcode';
                end
            end
            %Determine whether we are generating C or M files
            if(strcmpi(opts.cbmode,'mcode'))
                fgen = 'M';
            elseif(strcmpi(opts.cbmode,'cppad'))
                fgen = 'A';
            else
                fgen = 'C';
            end
            %Create New Callback Function based on req'd callbacks
            opts.havCon = any(B.conlin > 0);
            if(fgen == 'M')
                B.buildMFun('new',[],[],opts);
            else
                B.buildCFun('new',[],[],opts);
            end
            %Get Objective (Nonlinear or not, returned as f(x))
            ind = B.objlin >= 1; 
            if(sum(ind) == 0)
                error('This object does not contain an objective');
            end
            %Build Objective Callback
            if(fgen == 'M')
                prob.fun = B.buildMFun('obj',B.sobj(ind),B.vars,opts);
            else
                prob.fun = B.buildCFun('obj',B.sobj(ind),B.vars,opts);
            end
            %If requested, build gradient callback
            if(strcmpi(opts.use1stDerivs,'yes'))
                if(fgen == 'A')
                    prob.f = @(x) symb_ccb('grad',x); %no need to generate
                else
                    ajac = B.jac(ind,:);
                    if(fgen == 'M')
                        prob.f = B.buildMFun('grad',ajac,B.vars,opts);
                    else
                        prob.f = B.buildCFun('grad',ajac,B.vars,opts);
                    end
                end
            elseif(fgen=='C')
                B.buildCFun('grad'); %default empty
            end

            %If we have constraints, create callbacks
            if(opts.havCon)
                %Get ALL Constraint Indices
                ind = B.conlin > 0;                
                %Build Con Callback + row constraints
                if(fgen == 'M')
                    prob.nlcon = B.buildMFun('con',B.sobj(ind),B.vars,opts);
                else
                    prob.nlcon = B.buildCFun('con',B.sobj(ind),B.vars,opts);
                end
                prob.cl = B.cl(ind);
                prob.cu = B.cu(ind);
                %If we want derivatives, create jacobian + structure
                if(strcmpi(opts.use1stDerivs,'yes'))
                    if(fgen == 'A')
                        prob.nljac = @(x) symb_ccb('jac',x); %no need to generate
                        prob.nljacstr = @() symb_ccb('jacstr');
                    else
                        ajac = B.jac(ind,:);
                        if(fgen == 'M')
                            prob.nljac = B.buildMFun('jac',ajac,B.vars,opts);
                        else
                            prob.nljac = B.buildCFun('jac',ajac,B.vars,opts);
                        end  
                        %Build Jacobian Structure Callback
                        jacstr = zeros(size(ajac));
                        jacstr(logical(ajac ~= 0)) = 1;
                        prob.nljacstr = @() sparse(jacstr);
                    end                    
                elseif(fgen=='C')
                    B.buildCFun('jac'); %default empty
                end
            elseif(fgen=='C' || fgen=='A')
                B.buildCFun('con',[],[],opts); %default empty
                B.buildCFun('jac',[],[],opts); %default empty
            end
            %If we want 2nd derivatives, add them too
            if(strcmpi(opts.use2ndDerivs,'yes'))
                if(fgen == 'A')
                    prob.H = @(x,sigma,lambda) symb_ccb('hess',x,sigma,lambda); %no need to generate
                    prob.Hstr = @() symb_ccb('hstr');
                else
                    %Check if we have done a full build
                    if(isempty(B.hess))
                        if(B.verbose), optiwarn('No 2nd Derivatives Added - You may only obtain the Hessian of the Lagrangian if you called Build().\nDraft() only computes 1st Derivatives%s\n','.'); end
                        if(fgen=='C')
                            B.buildCFun('hess'); %default empty
                        end
                    else
                        %Build HessLag
                        buildHessLag(B);
                        %Build NL Hessian Callback
                        if(fgen=='M')
                            prob.H = B.buildMFun('hess',B.hesslag,B.vars,opts);
                        else
                            prob.H = B.buildCFun('hess',B.hesslag,B.vars,opts,B.noCons);
                        end
                        %Build Hessian Structure
                        lHstr = zeros(size(B.hesslag));
                        lHstr(logical(B.hesslag ~= 0)) = 1;
                        prob.Hstr = @() sparse(lHstr);
                    end
                end
            elseif(fgen=='C')
                B.buildCFun('hess'); %default empty
            end 

            %Bounds & Integer Vars
            prob.lb = B.lb;
            prob.ub = B.ub;
            prob.int = B.xtype; 
            %Problem Name
            if(any(B.xtype=='I'))
                prob.Name = 'SymBuilder MINLP';
            else
                prob.Name = 'SymBuilder NLP';
            end       
            
            %Generate pretend x0 incase none given
            prob.x0 = (prob.ub - prob.lb)./2;
            prob.x0(isnan(prob.x0)) = 0;
            prob.x0(isinf(prob.x0)) = 0;
            
            %Compile if C Code
            if(fgen ~= 'M')
                if(B.verbose), fprintf('Compiling....'); end
                %If compiling with VS2015 but pre R2015b, need to manually add in UCRT location
                cc = mex.getCompilerConfigurations(); post = [];
                for i = 1:length(cc)
                    if(~isempty(strfind(cc(i).Name,'Microsoft Visual C++')) && str2double(cc(i).Version) >= 14 && verLessThan('matlab','8.6'))
                        post = opti_FindUCRT();
                        break;
                    end
                end
                if(fgen == 'A')
                    src = 'symb_ccb.cpp';
                    %find cppad source (expected in utilities folder)
                    str = which('asl');
                    if(isempty(str)), error('Cannot find asl to locate CppAD source!'); end
                    idx = strfind(str,filesep); str = str(1:idx(end));
                    str = [str 'Source'];
                    str = ['mex -largeArrayDims symb_ccb.cpp -I"' str '" ' post];
                    if(strcmpi(computer,'PCWIN'))
                        str = [str ' -DCPPAD_SIZE_T_SAME_UNSIGNED_INT=1'];
                    end
                	eval(str);
                else
                    src = 'symb_ccb.c';
                    str = ['mex -largeArrayDims symb_ccb.c ' post];
                    eval(str);
                end
                if(B.verbose), fprintf('Done\n'); end
                if(strcmpi(opts.srckeep,'no'))
                    delete(src);
                end
            end
        end          
                
        %Get Build Status
        function tf = IsBuilt(B,doErr)
            if(nargin < 2), doErr = false; end               
            if(strcmp(B.bldstat,'unbuilt'))
                tf = false;
                if(doErr)
                    error('Please build the object using Draft() or Build()');
                end
            else
                tf = true;
            end
        end
        
        function [x,fval,ef,info] = solve(B,x0,opts)
        %Lowercase version
            if(nargin < 3), opts = []; end
            if(nargin < 2), x0 = []; end
            [x,fval,ef,info] = Solve(B,x0,opts);
        end
        
        function [x,fval,ef,info] = Solve(B,x0,opts)
        %Solve Object using OPTI

            %Check is built
            if(~IsBuilt(B))
                error('Please build the object first using Draft() or Build()');
            end
            if(nargin < 3 || isempty(opts)), opts = symbset; end
            if(nargin < 2), x0 = []; end
            
            if(B.verbose), fprintf('\nGenerating OPTI Object....\n'); end
            %Cusomtize settings based on solver
            switch(lower(opts.solver))
                %White box solvers
                case {'scip','baron'}
                     opts = symbset(opts,'use1stDerivs','no','use2ndDerivs','no','preallocate','no');
                %Derivative free solvers
                case {'nomad','pswarm','gmatlab'}                    
                    opts = symbset(opts,'use1stDerivs','no','use2ndDerivs','no');
                %No Hessian Support                   
                case {'filtersd','lbfgsb','nlopt'}
                    opts = symbset(opts,'use2ndDerivs','no');
            end
            B.Opt = GetOPTI(B,opts);
            if(B.verbose), fprintf('Done\n\n'); end
            %Solve
            [x,fval,ef,info] = solve(B.Opt,x0);
        end
        
        function B = plus(B,C)
        %Add objective, constraint, bound or integer declaration to SymBuilder Object
            if(~isa(B,'SymBuilder'))
                error('The LHS object must be a SymBuilder Object');
            end
            if(~ischar(C))
                error('The RHS object must be a string');
            end
            %Check for objective
            idx = strfind(C,'min');
            if(~isempty(idx))
                if(length(idx) > 1), error('min most only appear once in the objective'); end
                B.AddObj(C(idx(1)+4:end));
                return;
            end
            idx = strfind(C,'max');
            if(~isempty(idx))
                if(length(idx) > 1), error('max most only appear once in the objective'); end
                B.AddObj(['-(' C(idx(1)+4:end) ')']);
                return;
            end
            %Check for bound or integer (no mathematical ops)
            idx = regexp(C,{'*','+','-','/','^'});
            if(~isempty([idx{:}]))               
                B.AddCon(C);
            else
                idx = regexp(C,{'<=','>='});
                if(~isempty([idx{:}]))
                    B.AddBound(C);
                else
                    B.AddInteger(C);
                end
            end
                
        end
    end
    

    methods(Access=private)
                    
        function buildSymRep(B)
        %Build symbolic representation of the equations
        
            %Check if something new has been added
            if(strcmp(B.bldstat,'unbuilt'))
                %Build Symbolic Object                
                symobj = B.eqs;
                %Substitute expressions
                if(~isempty(B.exprsn))                                           
                    %Now subs into full equation system
                    wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
                    symobj = subs(symobj,B.exprsn(:,1),B.exprsn(:,2));
                    warning(wstate);
                    %Have to repeat until all nested expressions are sub'd
                    n = 10;  %max depth
                    no = size(B.exprsn,1);
                    se = sym(B.exprsn(:,1));
                    while n > 0
                        v = symvar(symobj); alldone = 1;
                        for i = 1:no
                            if(any(se(i) == v))
                                wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
                                symobj = subs(symobj,B.exprsn(i,1),B.exprsn(i,2));
                                warning(wstate);
                                alldone = 0;
                            else
                                if(i == no && alldone) %ensure we have checked them all
                                    n = -10;
                                    break;
                                end
                            end
                        end  
                        n = n - 1;
                    end
                    if(n == 0)
                        error('Maximum expression recursion depth reached!');
                    else
%                         fprintf('%d exp recursions required\n',n+10);
                    end
                end
                %Convert constants to decimal form
%                 digits(16); %set vpa digits
%                 BV = cell(size(B.constnt,1),1);
%                 for i = 1:size(B.constnt,1)
%                     BV{i} = sym(B.constnt{i,2},'d');
%                 end
                %Substitute constants
                if(~isempty(B.constnt))
                    wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
                    symobj = subs(symobj,B.constnt(:,1),B.constnt(:,2));
                    warning(wstate);
                end
                %Save symbolic vector
                B.sobj = symobj;
                %Save Variables
                B.vars = symvar(B.sobj);
                
                %Build default bounds and xtype
                n = length(B.vars);
                B.lb = -1*Inf(n,1);
                B.ub = Inf(n,1);
                B.xtype = repmat('C',1,n);   
                %Get names of variables and their indices
                names = SymBuilder.detVarNames(B.vars);
                if(isempty(names) && size(B.bnds,1) > 0)
                    error('Cannot process bounds without declared variables!');
                end
                                
                %Process bound declarations
                for i = 1:size(B.bnds,1)
                    %Check if bounds are numerical
                    llb = str2double(B.bnds{i,2});
                    if(isnan(llb)) %not a number, check constants
                        try
                            ind = strcmp(B.constnt(:,1),B.bnds{i,2});
                            if(any(ind))
                                llb = B.constnt{ind,2};
                            else
                                error('Unknown lower bound: %s',B.bnds{i,2});
                            end
                        catch
                            error('Error processing lower bound for var %d',i);
                        end
                    end
                    lub = str2double(B.bnds{i,3});
                    if(isnan(lub)) %not a number, check constants
                        try
                            ind = strcmp(B.constnt(:,1),B.bnds{i,3});
                            if(any(ind))
                                lub = B.constnt{ind,2};
                            else
                                error('Unknown upper bound: %s',B.bnds{i,3});
                            end
                        catch
                            error('Error processing upper bound for var %d',i);
                        end
                    end
                    %Check if variable exists
                    ind = logical(B.vars == sym(B.bnds{i,1}));
                    if(any(ind))
                        if(~isinf(llb)), B.lb(ind) = llb; end
                        if(~isinf(lub)), B.ub(ind) = lub; end                    
                    else %could be applied to multiple variables
                        ind = strcmp(names(:,1),B.bnds{i,1});
                        if(isempty(ind) || all(ind  == 0))
                            error('Unknown bound: %s',B.bnds{i,1});
                        elseif(sum(ind) > 1)
                            error('Bound name is ambiguous and matches more than one variable name');
                        else                            
                            %Extract indices
                            start = names{ind,2}(1);
                            send = names{ind,2}(end);
                            if(~isinf(llb)), B.lb(start:send) = llb; end
                            if(~isinf(lub)), B.ub(start:send) = lub; end 
                        end
                    end
                end
                
                %Process integer declarations
                for i = 1:size(B.vartypes,1)
                    %Check if variable exists
                    ind = logical(B.vars == sym(B.vartypes{i,1}));
                    if(any(ind))
                         B.xtype(ind) = B.vartypes{i,2};                
                    else %could be applied to multiple variables
                        ind = strcmp(names(:,1),B.vartypes{i,1});
                        if(isempty(ind))
                            error('Unknown variable: %s',B.vartypes{i,1});
                        elseif(sum(ind) > 1)
                            error('Variable name is ambiguous and matches more than one variable name');
                        else                            
                            %Extract indices
                            start = names{ind,2}(1);
                            send = names{ind,2}(end);
                            B.xtype(start:send) = B.vartypes{i,2};
                        end
                    end
                end
            end
        end
                
        function buildJac(B)
        %Evaluate Jacobian of entire equation system
        
            if(isempty(B.sobj))
                error('You must build a symbolic representation of the system first!');
            end
            %Call SymToolbox jacobian function
            B.jac = jacobian(B.sobj);
        end  
                
        function buildHess(B)
        %Evaluate Hessian(s) of entire equation system
        
            if(isempty(B.jac))
                error('You must evaluate the Jacobian of the system first!');
            end
            noeq = B.noObjs + B.noCons;
            B.hess = struct('H',cell(noeq,1),'lin',cell(noeq,1));            
            %Build Index Vectors
            B.conlin = zeros(noeq,1);
            B.objlin = zeros(noeq,1);
            B.objbias = zeros(noeq,1);
            %Go through each equation
            for i = 1:noeq
                J = B.jac(i,:);
                %Determine if linear
                if(isempty(symvar(J)))
                    B.hess(i).lin = 1;
                    B.hess(i).H = 0;
                else
                    %Hessian is the Jacobian of the Jacobian (in this case)
                    H = jacobian(J,B.vars);
                    s = symvar(H);
                    %Determine if quadratic
                    if(isempty(s))
                        B.hess(i).lin = 2;
                        B.hess(i).H = H;
                    else
                        %Otherwise must be nonlinear, save entire symobj (for now)
                        B.hess(i).lin = 3;
                        B.hess(i).H = H;
                    end
                end
                %If linear or quadratic, save bias term
                if(B.hess(i).lin <= 2)
                    eq = B.sobj(i);
                    v = symvar(eq);
                    %Find bias term
                    B.objbias(i) = double(subs(eq,v,{zeros(size(v))}));
                end 
                    
                %Save linearity
                if(any(B.indobj == i))                       
                    B.objlin(i) = B.hess(i).lin;
                else
                    B.conlin(i) = B.hess(i).lin;
                end
            end
        end
                
        function buildHessLag(B)
        %Build Lagrangian Hessian assuming IPOPT (tril) format 
        
            if(isempty(B.hess))
                error('You must evaluate the Hessian of the system first!');
            elseif(~isempty(B.hesslag))
                return; %already done
            end
            %Start with objective
            nvars = length(B.vars);
            ind = B.objlin >= 2;
            if(any(ind)) %check we have a nonlinear (or quadratic) objective  
                if(sum(ind) > 1)
                    error('This interface only support single objective problems');
                end
                s = sym('sigma');
                objH = B.hess(ind).H;
                H = s*objH;
            else
                H = zeros(nvars);
            end
            %Now each nonlinear (or quadratic) constraint
            ind = B.conlin > 1;
            %Two sets of indexes, one for lambda, one for B.hess
            index_BH = find(ind);
            ind(B.indobj) = [];
            index_L = find(ind);
            wstate = warning('off','symbolic:sym:sym:DeprecateExpressions');
            for i = 1:length(index_BH)
                l = sym(sprintf('lambda(%d)',index_L(i)));
                H = H + l*B.hess(index_BH(i)).H;
            end                  
            warning(wstate);
            %Save resulting Hessian
            B.hesslag = tril(H);
        end
                
        function detLinearity(B)
        %Determine linearity of the supplied equations (without Hessian)
        
            if(isempty(B.jac))
                error('You must evaluate the Jacobian of the system first!');
            end
            %Build Index Vectors
            noeq = B.noCons + B.noObjs;
            B.conlin = zeros(noeq,1);
            B.objlin = zeros(noeq,1);
            %Determine Linearity (based on existence of symvars)
            for i = 1:noeq %must search all equations
                %Gather first derivative vars
                s = symvar(B.jac(i,:));
                %Objective
                if(any(B.indobj == i))
                    if(isempty(s))
                        B.objlin(i) = 1; %lin
                    else
                        B.objlin(i) = 3; %nonlin
                    end
                else                   
                    if(isempty(s))
                        B.conlin(i) = 1;
                    else
                        B.conlin(i) = 3;
                    end
                end
            end            
        end 
    end
    
    
    methods(Static)        
        %Parse Constraint String
        [str,l,u] = parseConstraint(str);       
        %Extract RHS of expression string, remove from eq, and return LHS
        [str,exp] = parseExpression(str);
        %Extract Group and Name from Result Expression Name
        [group,name] = parseResExp(name);
        %Parse bound string
        [var,lb,ub] = parseBound(str);
        %Parse integer string
        [var,xtype] = parseInteger(str);
        %Determine Variable names
        varn = detVarNames(svar);
        %Build MATLAB function file
        cb = buildMFun(mode,sobj,svar,opts);
        %Build C Code function file
        cb = buildCFun(mode,sobj,svar,opts,nocon);
        %Convert symbolic expression into matlab function
        fun = sym2fun(sobj,svar,var,skipSubs);     
        
        %Check we have a compiler suitable for use with Cppad
        function ok = CheckCppADCompile(doerr)
            ok = true;
            c = mex.getCompilerConfigurations('c++');
            if(isempty(c) || (isempty(strfind(c.Name,'Microsoft Visual C++')) && isempty(strfind(c.Name,'Intel C++'))))
                str = sprintf(['CppAD compilation is only available with Microsoft Visual Studio or Intel as the C++ compiler.\n\nPlease either select it'...
                           ' as the MEX compiler (via mex -setup), or download from Microsoft:\n\n %s (2012, recommended)\n\nor\n\n %s (2013)'],...
                           'http://www.microsoft.com/en-nz/download/details.aspx?id=34673','http://www.visualstudio.com/en-us/products/visual-studio-express-vs.aspx');
                if(doerr)
                    error(str) %#ok<SPERR>
                else %warning, auto has selected it
                    ok = false;
                    optiwarn('symb:cppad',sprintf('SymBuilder has chosen CppAD has the most effective callback strategy, but %s',str));
                end
            end
        end
    end
    
end

