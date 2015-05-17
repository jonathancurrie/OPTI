function displayOPTI(O)
%Display OPTI Problem Information
%
%   Called By OPTI Class

%   Copyright (C) 2011 Jonathan Currie (I2C2)

%Distribute Structure Variables;
prob = O.prob;
opts = O.opts;
sizes = prob.sizes;
nineq = sizes.nineq;
neq = sizes.neq;
nbnds = sizes.nbnds;
nqc = sizes.nqc;
nsdcone = sizes.nsdcone;
nnlineq = sizes.nnlineq;
nnleq = sizes.nnleq;
nsos1 = sizes.nsos1;
nsos2 = sizes.nsos2;
sense = prob.sense;

%Types of constraints
lin = 'rl <= Ax <= ru';
linineq = 'Ax <= b';
lineq = 'Aeqx = beq';
bnds = 'lb <= x <= ub';
quad = 'qrl <= x''Qx + l''x <= qru';
sdp = 'F1*x1 + F2*x2 + ... + Fn*xn - F0 >= 0 [PSD]';
sos1 = 'x in SOS1';
sos2 = 'x in SOS2';
nlineq = 'c(x) <= d';
nleq = 'ceq(x) = deq';
nl = 'cl <= c(x) <= cu';
str1 = ' s.t. '; 
str2 = '      '; 

%Constant objective term
if(isfield(prob,'objbias') && ~isempty(prob.objbias) && prob.objbias ~= 0)
    objc = '+ objbias';
else
    objc = '';
end

%See if we are solving as something other than problem definition
if(~strcmpi(prob.type,prob.solvetype))
    stype = sprintf(' (Solved as %s)',upper(prob.solvetype));
else
    stype = '';
end

%**************** Title ****************%
disp('------------------------------------------------------');
switch(prob.type)
    case 'SLE'
        disp('System of Linear Equations (SLE)');
        disp(' Ax = b');
    case 'LP'
        disp('Linear Program (LP) Optimization');
        if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end            
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end
    case 'MILP'
        disp('Mixed Integer Linear Program (MILP) Optimization');
        if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end        
        if(sizes.nint || sizes.nbin)
            disp('      xi = Integer / Binary');
        end
        if(nsos1),   fprintf('%s%s\n',str1,sos1); end
        if(nsos2),   fprintf('%s%s\n',str1,sos2); end
    case 'BILP'
        disp('Binary Integer Linear Program (BILP) Optimization');
        if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end
        disp('      xb = Binary');
    case 'QP'
        if(prob.iscon)
            fprintf('Quadratic Program (QP) Optimization%s\n',stype);
            if(sense==1); fprintf(' min  0.5x''Hx + f''x %s\n',objc); else fprintf(' max  1/2x''Hx + f''x %s\n',objc); end            
            if(~isempty(prob.rl))
                fprintf('%s%s\n',str1,lin); str1 = str2;
            else
                if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
            end
            if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
            if(nbnds),  fprintf('%s%s\n',str1,bnds); end
        else
            disp('Unconstrained Quadratic Program (UQP)');
            if(sense==1); disp(' min  1/2x''Hx + f''x'); else disp(' max  1/2x''Hx + f''x'); end
        end
    case 'QCQP'
        if(nnz(prob.H)==0)
            fprintf('Quadratically Constrained Linear Program (QCLP) Optimization%s\n',stype);        
            if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        else
            fprintf('Quadratically Constrained Quadratic Program (QCQP) Optimization%s\n',stype);        
            if(sense==1); fprintf(' min  0.5x''Hx + f''x %s\n',objc); else fprintf(' max  1/2x''Hx + f''x %s\n',objc); end 
        end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nqc),    fprintf('%s%s\n',str1,quad); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end
    case 'MIQP'
        fprintf('Mixed Integer Quadratic Program (MIQP) Optimization%s\n',stype);
        if(sense==1); fprintf(' min  0.5x''Hx + f''x %s\n',objc); else fprintf(' max  1/2x''Hx + f''x %s\n',objc); end 
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end    
        disp('      xi = Integer / Binary');
    case 'MIQCQP'
        if(nnz(prob.H)==0)
            fprintf('Mixed Integer Quadratically Constrained Linear Program (MIQCLP) Optimization%s\n',stype);
            if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end 
        else
            fprintf('Mixed Integer Quadratically Constrained Quadratic Program (MIQCQP) Optimization%s\n',stype);
            if(sense==1); fprintf(' min  0.5x''Hx + f''x %s\n',objc); else fprintf(' max  1/2x''Hx + f''x %s\n',objc); end 
        end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nqc),    fprintf('%s%s\n',str1,quad); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end    
        disp('      xi = Integer / Binary');
    case 'SDP'
        disp('Semidefinite Program (SDP) Optimization');
        if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end            
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nsdcone),    fprintf('%s%s\n',str1,sdp); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end
    case 'MISDP'
        disp('Mixed Integer Semidefinite Program (MISDP) Optimization');
        if(sense==1); fprintf(' min  f''x %s\n',objc); else fprintf(' max  f''x %s\n',objc); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end            
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nsdcone),    fprintf('%s%s\n',str1,sdp); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end       
        disp('      xi = Integer / Binary');
    case 'UNO'
        disp('Unconstrained Nonlinear Optimization (UNO)');
        if(sense==1); disp(' min  f(x)'); else disp(' max  f(x)'); end
    case 'SNLE'
        fprintf('System of Nonlinear Equations (SNLE)%s\n',stype);
        disp(' F(x) = 0');
    case 'SCNLE'
        fprintf('System of Constrained Nonlinear Equations (SCNLE)%s\n',stype);
        disp(' F(x) = 0');
        if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end 
        if(sizes.nint || sizes.nbin)
            disp('      xi = Integer / Binary');
        end
    case 'NLS'
        if(~isempty(prob.weighting))
            fprintf('Weighted Nonlinear Least Squares (NLS)%s\n',stype);
            wts = true;
        else
            fprintf('Nonlinear Least Squares (NLS)%s\n',stype);
            wts = false;
        end
        if(~isempty(prob.xdata))
            if(wts)
                s = '(F(x,xdata) - ydata).*sqrt(wts)';
            else
                s = 'F(x,xdata) - ydata';
            end
        else
            if(wts)
                s = '(F(x) - ydata).*wts';
            else
                s = 'F(x) - ydata';
            end
        end            
        if(sense==1)
            fprintf(' min sum[( %s ).^2]\n',s);
        else
            fprintf(' max sum[( %s ).^2]\n',s);
        end
        if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end      
    case 'DNLS'
        if(~isempty(prob.weighting))
            fprintf('Weighted Dynamic Parameter Estimation Problem (DNLS)%s\n',stype);
            wts = true;
        else
            fprintf('Dynamic Parameter Estimation Problem (DNLS)%s\n',stype);
            wts = false;
        end        
        if(~isempty(prob.xdata))
            if(wts)
                disp(' min sum[( (int[F(t,z,p)] - ydata).*sqrt(wts) ).^2]');
            else
                disp(' min sum[( int[F(t,z,p)] - ydata ).^2]');
            end
        end
        if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); end 
    case 'NLP'
        disp('Nonlinear Program (NLP) Optimization');
        if(sense==1); disp(' min  f(x)'); else disp(' max  f(x)'); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
            if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); str1 = str2; end
        if(~isempty(prob.nlrhs))
            if(nnlineq),fprintf('%s%s\n',str1,nlineq); str1 = str2; end
            if(nnleq),  fprintf('%s%s\n',str1,nleq); end
        elseif(~isempty(prob.cl))
            fprintf('%s%s\n',str1,nl);
        end
    case 'MINLP'
        disp('Mixed Integer Nonlinear Program (MINLP) Optimization');
        if(sense==1); disp(' min  f(x)'); else disp(' max  f(x)'); end
        if(~isempty(prob.rl))
            fprintf('%s%s\n',str1,lin); str1 = str2;
        else
           if(nineq),  fprintf('%s%s\n',str1,linineq); str1 = str2; end
        end
        if(neq && ~isempty(prob.beq)),    fprintf('%s%s\n',str1,lineq); str1 = str2; end
        if(nbnds),  fprintf('%s%s\n',str1,bnds); str1 = str2; end
        if(~isempty(prob.nlrhs))
            if(nnlineq),fprintf('%s%s\n',str1,nlineq); str1 = str2; end
            if(nnleq),  fprintf('%s%s\n',str1,nleq); end
        elseif(~isempty(prob.cl))
            fprintf('%s%s\n',str1,nl);
        end
        disp('      xi = Integer / Binary');
end

%**************** Properties ****************%
disp('------------------------------------------------------');
disp('   Problem Properties: ')
if(isempty(sizes.ndec))
    fprintf('# Decision Variables:     Cannot Determine\n');   
else
    if(~isempty(prob.H) && issparse(prob.H) && nnz(prob.H)>0)
        fprintf('# Decision Variables:      %3d [H: %d nz]\n',sizes.ndec,nnz(prob.misc.Horig)); %QP Hessian
    else
        fprintf('# Decision Variables:      %3d\n',sizes.ndec);
        if(strcmpi(prob.type,'DNLS'))
            fprintf('  # Parameters:            %3d\n',length(prob.x0)-sum(isnan(prob.odez0)));
            fprintf('  # Initial Conditions:    %3d\n',sum(isnan(prob.odez0)));
        end
    end
end
if(~isempty(prob.ydata))
    fprintf('# Data Points:             %3d\n',numel(prob.ydata));
end
if(strcmpi(prob.type,'SLE'))
    fprintf('               LHS A:      %d x %d [A: %d nz]\n',size(prob.A,1),size(prob.A,2),nnz(prob.A));
    fprintf('               RHS b:      %d x 1\n',length(prob.b));
elseif(prob.iscon)
    fprintf('# Constraints:             %3d\n',sizes.ncon);
    if(nineq)
        if(~issparse(prob.A) || ~isempty(prob.rl))
            fprintf('  # Linear Inequality:     %3d\n',sizes.nineq);
        else
            fprintf('  # Linear Inequality:     %3d [A: %d nz]\n',sizes.nineq,nnz(prob.A));
        end
    end
    if(neq)
        if(~issparse(prob.Aeq))
            fprintf('  # Linear Equality:       %3d\n',sizes.neq);
        else
            fprintf('  # Linear Equality:       %3d [Aeq: %d nz]\n',sizes.neq,nnz(prob.Aeq));
        end
    end
    if(nqc)
        if(~issparse(prob.Q))
            fprintf('  # Quadratic Constraints: %3d\n',sizes.nqc);
        else 
            fprintf('  # Quadratic Constraints: %3d [Q: %d nz]\n',sizes.nqc,nnz(prob.Q));
        end
    end
    if(nsdcone)
        if(~issparse(prob.sdcone))
            fprintf('  # Semidefinite Cones:    %3d\n',sizes.nsdcone);
        else 
            fprintf('  # Semidefinite Cones:    %3d [SDCone: %d nz]\n',sizes.nsdcone,nnz(prob.sdcone));
        end
    end
    if(sizes.nbnds)
        fprintf('  # Bounds:                %3d\n',sizes.nbnds);
    end
    
    if(sizes.nbin)
        fprintf('  # Binary Variables:      %3d\n',sizes.nbin);
    end
    if(sizes.nint)
        fprintf('  # Integer Variables:     %3d\n',sizes.nint);
    end
    if(nsos1)
        fprintf('  # SOS1:                  %3d\n',sizes.nsos1);
    end
    if(nsos2)
        fprintf('  # SOS2:                  %3d\n',sizes.nsos2);
    end
    
    % IF SPARSE

    if(sizes.nnlineq)
        fprintf('  # Nonlinear Inequality:  %3d\n',sizes.nnlineq);
    end
    if(sizes.nnleq)
        fprintf('  # Nonlinear Equality:    %3d\n',sizes.nnleq);
    end
end

disp('------------------------------------------------------');

%**************** Solver ****************%
disp('  Solver Parameters:')
if(prob.ampl.useASL)
    fprintf('Solver:                    %s [Using AMPL Interface]\n',upper(opts.solver));    
else
    dopts = optidynset(opts.dynamicOpts);
    if(strcmpi(opts.solver,'nlopt'))
        fprintf('Solver:                    %s with %s\n',upper(opts.solver),upper(nloptSolver(O.nlprob.options.algorithm)));
    else
        fprintf('Solver:                    %s\n',upper(opts.solver));
    end
    
    if(strcmpi(prob.type,'DNLS'))
        fprintf('ODE Integrator:            %s\n',upper(dopts.integrator));        
        if(~strcmpi(dopts.sensitivity,'none')) 
            if(isempty(dopts.sensitivity)), dopts.sensitivity = 'user'; end
            
            switch(lower(dopts.sensitivity))
                case 'nd', dfdz = 'Numerical Differentiation';
                case 'ad', dfdz = 'Automatic Differentiation';   
                case 'cs', dfdz = 'Complex-Step Differentiation';
                case 'user'
                    if(~isempty(dopts.dfdz))
                        dfdz = 'User Supplied';
                    else
                        dfdz = 'Numerical Differentiation'; %default alternative
                    end
            end            
            
            switch(lower(dopts.sensitivity))
                case 'nd', dfdp = 'Numerical Differentiation';
                case 'ad', dfdp = 'Automatic Differentiation';    
                case 'cs', dfdp = 'Complex-Step Differentiation';
                case 'user'
                    if(~isempty(dopts.dfdz))
                        dfdp = 'User Supplied';
                    else
                        dfdp = 'Numerical Differentiation'; %default alternative
                    end
            end 
            fprintf('Sensitivity DFDZ:          %s\n',dfdz);
            fprintf('Sensitivity DFDP:          %s\n',dfdp);
        end

    elseif(any(strcmpi(prob.type,{'SNLE','SCNLE','NLS','NLP','UNO','MINLP'})) && isa(prob.fun,'function_handle'))
        if(prob.numdif.grad), fnd = ' [numdiff]'; else fnd = ''; end
        %Get solver info
        info = optiSolverInfo(opts.solver,prob.type,O.nlprob,opts);
        if(info.der1 || (strcmpi(opts.solver,'matlab') && strcmpi(prob.type,'NLP'))) %matlab is optional
            derstr = 'Not Supplied'; 
        else
            derstr = 'Not Required'; 
        end
        if(isempty(prob.f))
            f = derstr;
        elseif(length(func2str(prob.f)) > 60)
            f = ['@(x) ...' fnd];
        else
            f = [func2str(prob.f) fnd];
        end              
        if(prob.numdif.hess), Hnd = ' [numdiff]'; else Hnd = ''; end       
        %Check if solver can use second derivs & not vector function
        if(info.der2 && ~any(strcmpi(prob.type,{'snle','nls'}))) 
            if(isempty(prob.H))
                H = derstr;
            elseif(length(func2str(prob.H)) > 60)
                H = ['@(x) ...' Hnd];
            else
                H = [func2str(prob.H) Hnd];
            end
            if(isempty(prob.Hstr))
                Hs = derstr;
            else
                Hs = 'Supplied';
            end
        else
            H = []; Hs = [];
        end
        if(~isempty(prob.nlcon))
            if(prob.numdif.jac), jnd = ' [numdiff]'; else jnd = ''; end
            if(isempty(prob.nljac))
                j = derstr;
            elseif(length(func2str(prob.nljac)) > 60)
                j = ['@(x) ...' jnd];
            else
                j = [func2str(prob.nljac) jnd];
            end
            if(isempty(prob.nljacstr))
                js = derstr;
            else
                js = 'Supplied';
            end
        else
            j = []; js = [];
        end
        
        if(prob.misc.forceSNLE && ~strcmpi(prob.solvetype,'nlp') && ~strcmpi(prob.solvetype,'minlp'))
            fprintf('Equation Jacobian:         %s\n',f);
        else
            fprintf('Objective Gradient:        %s\n',f);
            if(~isempty(j)), fprintf('Constraint Jacobian:       %s\n',j); end
        end        
        if(~isempty(js))
            if(strcmpi(js,'Not Supplied'))
                fprintf('Jacobian Structure:        Assumed 100%% Dense\n'); 
            else
                fprintf('Jacobian Structure:        %s\n',js); 
            end
        end
        if(~isempty(H))
            if(isempty(j))
                fprintf('Objective Hessian:         %s\n',H); 
            else
                fprintf('Lagrangian Hessian:        %s\n',H);
            end
        end
        if(~isempty(Hs)), fprintf('Hessian Structure:         %s\n',Hs); end
    end
end

disp('------------------------------------------------------'); 
  
    
end