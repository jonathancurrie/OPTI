function [x,fval,exitflag,info] = opti_scipnl(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xint,x0,opts)
%OPTI_SCIPNL Solve a NLP/MINLP using SCIP to Global Optimality
%
%   min fun(x)                 subject to:     rl <= A*x <= ru
%    x                                         lb <= x <= ub
%                                              cl <= nlcon(x) <= cu
%                                              for i = 1..n: xi in Z
%                                              for j = 1..m: xj in {0,1} [i!=j]
%
%   Full Calling Form:
%     [x,fval,exitflag,info] = opti_scipnl(fun,A,rl,ru,lb,ub,nlcon,cl,cu,xint,x0,opts)
%
%   NOTE:
%     This solver (and interface) parses your supplied function(s) into an
%     algebraic description of the problem. This means your functions must be
%     deterministic (no random numbers or stochastic elements), and must also
%     contain only a subset of functions (log,log10,exp,abs,sqrt,norm). In
%     addition, only MATLAB code may be parsed, meaning no Simulink or MEX
%     models may be optimized using this solver.
%
%
%   x = opti_scipnl(fun) solves an unconstrained NLP (UNO) with objective
%   specified by fun. Note you must supply x0 if solving an UNO.
%
%   x = opti_scipnl(fun,A,rl,ru) solves a NLP subject to linear constraints
%   where A is the linear constraint matrix, and rl and ru are lower and
%   upper row bounds.
%
%   x = opti_scipnl(fun,A,rl,ru,lb,ub) solves subject to upper and lower
%   decision variable bounds. Finite bounds on ALL decision variables are
%   HIGHLY recommended when solving global optimization problems.
%
%   x = opti_scipnl(fun,...,ub,nlcon,cl,cu) solves subject to nonlinear
%   constraints specified by nlcon (column vector of constraint functions)
%   and lower and upper row bounds (cl, cu).
%
%   x = opti_scipnl(fun,...,cu,xint) solves a MINLP subject to integer and
%   binary constraints specified by xint. Integer constraints are specified
%   using a string of integer variables ('C', 'I', 'B').
%
%   x = opti_scipnl(fun,...,xint,x0) uses x0 to supply a validation vector
%   of the decision variables for the MEX interface, to use to verify your 
%   objective and constraints have been successfully converted to SCIP
%   expressions. Note it DOES NOT supply a starting guess to the solver (currently).
%
%   x = opti_scipnl(fun,...,x0,opts) uses opts to pass optiset options to the
%   solver.
%
%   [x,fval,exitflag,info] = opti_scipnl(...) returns the objective value at
%   the solution, together with the solver exitflag, and an information
%   structure.
%
%   THIS IS A WRAPPER FOR SCIP USING THE MEX INTERFACE
%   See supplied ZIB Academic License

%   Copyright (C) 2012/2013 Jonathan Currie (I2C2)

t = tic;

% Handle missing arguments
if nargin < 12, opts = optiset; end
if nargin < 11, x0 = []; end
if nargin < 10, xint = []; end
if nargin < 9, cu = []; end
if nargin < 8, cl = []; end
if nargin < 7, nlcon = []; end
if nargin < 6, ub = []; end
if nargin < 5, lb = []; end
if nargin < 4, ru = []; end
if nargin < 3, rl = []; end
if nargin < 2, A = []; end
if nargin < 1, error('You must supply at least one argument to opti_scipnl'); end

%Check for number of decision variables
if(~isempty(x0)) %vector or matrix
    ndec = numel(x0);
elseif(~isempty(lb)) %vector
    ndec = length(lb);
elseif(~isempty(ub)) %vector
    ndec = length(ub);
elseif(~isempty(xint)) %vector
    ndec = length(xint);
elseif(~isempty(A))
    ndec = size(A,2);
else
    error('The MATLAB - SCIP Interface requires the number of variables to be specified.\n%s',...
          ['The interface cannot determine this number from your problem thus you will need',...
          ' to specify this via x0, e.g. x0 = ones(no_vars,1).']);
end
%Build x0 if empty, used for validation
if(isempty(x0))
    x0 = randn(ndec,1); 
end

warn = strcmpi(opts.warnings,'all');
%Check for NaNs
if(any(isnan(x0)))
    if(warn)
        optiwarn('opti:nan','For the SCIP interface x0 must not contain NaN. Replacing with random numbers.');
    end
    x0 = randn(ndec,1);
end

%Encode Nonlinear Constraints into SCIP MEX Interface instruction lists
x = scipvar(size(x0));
if(~isempty(nlcon))
    try
        n = nlcon(x);
    catch ME
        ex = processEqErr(ME,'a constraint'); 
        throwAsCaller(ex);
    end
    if(numel(n) > 1) %multiple constraints
        nl.instr = cell(numel(n),1);
        for i = 1:numel(n)            
            if(isnumeric(n(i))) %read as number
                nl.instr{i} = [0 n(i)]; 
            elseif(isempty(n(i).ins)) %assume just a variable
                nl.instr{i} = [1; n(i).indx];
            else
                nl.instr{i} = n(i).ins;
            end
        end        
    else        
        if(isnumeric(n)) %read as number
            nl.instr = [0 n]; 
        elseif(isempty(n.ins)) %assume just a variable
            nl.instr = [1; n.indx];
        else
            nl.instr = n.ins;
        end
    end
    nl.cl = cl;
    nl.cu = cu;
    %Verification Fields
    nl.nlcon_val = nlcon(x0);
    nl.x0 = x0;
    %Check for Inf or NaN
    if(any(isnan(nl.nlcon_val)) || any(isinf(nl.nlcon_val)))
        error('One or more constraints resulted in Inf or NaN at the initial guess (x0). Please provide a better initial guess vector.');
    end
end
%Encode Nonlinear Objective into SCIP MEX Interface instruction list
if(~isempty(fun))
    try
        f = fun(x);
    catch ME
        ex = processEqErr(ME,'the objective');
        throwAsCaller(ex);
    end    
    if(isnumeric(f)) %read as number
        nl.obj_instr = [0 f]; 
    elseif(isempty(f.ins)) %assume just a variable
        nl.obj_instr = [1; f.indx];
    else
        nl.obj_instr = f.ins;
    end
    %Verification field
    nl.obj_val = fun(x0);
    nl.x0 = x0;
    %Check for Inf or NaN
    if(any(isnan(nl.obj_val)) || any(isinf(nl.obj_val)))
        error('The objective resulted in Inf or NaN at the initial guess (x0). Please provide a better initial guess vector.');
    end
end
%Check sparsity
if(~isempty(A) && ~issparse(A))
    if(warn)
        optiwarn('opti:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    A = sparse(A);
end

%Setup Printing
opts.display = dispLevel(opts.display);

%Run SCIP
[x,fval,exitflag,stats] = scip([],zeros(ndec,1),A,rl,ru,lb,ub,xint,[],[],nl,opts);

%Reshape output
x = reshape(x,size(x0));

%Assign Outputs
info.BBNodes = stats.BBnodes;
info.BBGap = stats.BBgap;
info.Time = toc(t);
info.Algorithm = 'SCIP: Spatial Branch and Bound using IPOPT and SoPlex';

switch(exitflag)
    case 10
        if(any(xint ~= 'C'))
            info.Status = 'Globally Integer Optimal';
        else
            info.Status = 'Globally Optimal';
        end
        exitflag = 1;
    case {2,3,4,5,6,7,8}
        info.Status = 'Exceeded Iterations / Time / Nodes';
        exitflag = 0;
    case 11
        info.Status = 'Infeasible';
        exitflag = -1;
    case {12,13}
        info.Status = 'Unbounded or Infeasible';
        exitflag = -2;
    case 1
        info.Status = 'User Exited';
        exitflag = -5;
    otherwise
        info.Status = [];
end

function  print_level = dispLevel(lev)
%Return CLP compatible display level
switch(lower(lev))
    case'off'
        print_level = 0;
    case 'iter'
        print_level = 4;
    case 'final'
        print_level = 3;
end

function ex = processEqErr(ME,name)

str = [];
for i = 1:length(ME.stack)
    if(i==1)
        if(length(ME.stack)>1 && strcmpi(ME.stack(2).name,'opti_scipnl'))
            str = sprintf('%sError using "%s" (line %d)\n%s\n',str,ME.stack(i).name,ME.stack(i).line,ME.message);
        else
            str = sprintf('%sError using "%s" (line %d)\n%s\n\nError Trace:\n',str,ME.stack(i).name,ME.stack(i).line,ME.message);
        end
    else
        if(strcmpi(ME.stack(i).name,'opti_scipnl'))
            break;
        end
        str = sprintf('%s- Error in "%s" (line %d)\n',str,ME.stack(i).name,ME.stack(i).line);
    end
end

ex = MException('OPTI:SCIPVAR',['There was an error processing %s function into SCIP compatible form.\n'...
                                'Please examine the below error to correct your function:\n\n%s\nRemember '...
                                'SCIP only supports a subset of MATLAB commands. Use:\n>> methods(scipvar)\n\nto see compatible functions and operators.\n\n'],name,str);

                           