function [Opt,SymB] = optisym(fun,x0,lb,ub,con,cl,cu,xtype,sopts,verbose)
%OPTSYM  Generate an OPTI Optimization Problem using SymBuilder + Symbolic Toolbox
%
%  [optiObj,SymBobj] = optisym(fun,x0,lb,ub,con,cl,cu,xtype,sopts,verbose)
%
%   optiObj = optisym(fun,x0) takes a standard MATLAB function (fun) and
%   initial solution guess (x0) and converts the problem to a Symbolic
%   form. It then identifies the type of problem entered (L/Q/NLP) and
%   builds the corresponding OPTI object.
%
%   optiObj = optisym(fun,x0,lb,ub) allows the user to specify lower (lb)
%   and upper (ub) bounds on the problem.
%
%   optiObj = optisym(fun,x0,lb,ub,nlcon,cl,cu) allows the user to specify
%   a constraint function (con), together with lower (cl) and upper (cu)
%   constraint bounds. As with the objective function, the constraint
%   function will be automatically converted to Symbolic form and each
%   constraint identified as linear, quadratic or nonlinear.
%
%   optiObj = optisym(fun,...,cu,xtype) allows the user to specify variable
%   integrality via xtype, with 'C' as continuous, 'B' as binary and 'I' as
%   integer.
%
%   optiObj = optisym(fun,...,xtype,sopts) allows the user to specify
%   symbset options.
%
%   optiObj = optisym(fun,...,sopts,verbose) controls the verbosity of the
%   internal SymBuilder object. verbose = false (default) removes all
%   screen output.
%
%   [optiObj,SymBobj] = optisym(...) also returns the internal SymBuilder
%   object.
%
%   See also opti opti.solve 
%
%   Copyright (C) 2011-2014 Jonathan Currie (www.i2c2.aut.ac.nz)

if(nargin < 10), verbose = false; end
if(nargin < 9 || isempty(sopts)), sopts = symbset; end
if(nargin < 8), xtype = []; end
if(nargin < 7), cu = []; end
if(nargin < 6), cl = []; end
if(nargin == 5), error('You must supply nlcon, cl and cu for nonlinear constraints'); end
if(nargin < 5), con = []; end
if(nargin < 4), ub = []; end
if(nargin < 3), lb = []; end
if(nargin < 2), error('You must supply at least 2 arguments to this function'); end
if(isempty(which('syms'))), error('The Symbolic Toolbox must be installed to use this function'); end

if(isempty(x0))
    if(~isempty(lb))
        ndec = length(lb);
    elseif(~isempty(ub))
        ndec = length(ub);
    elseif(~isempty(xtype))
        ndec = length(xtype);
    else
        error('Please supply x0 in order to tell OPTI how many variables there are');
    end
    x0 = zeros(ndec,1); %just for error checking
else
    ndec = length(x0);
end

%Generate Decision Variable Vector
x = sym('x',[ndec,1]);

%Build SymBuilder Object
B = SymBuilder(verbose);

%Add Objective
try
    if(isa(fun,'sym'))
        B.AddObj(fun);
    else
        B.AddObj(fun(x));
    end
catch ME
    error('There was an error processing your objective into Symbolic Form. Please examine the error below and correct your function:\n\nError: %s',ME.message);
end

%Add Constraints
if(~isempty(con))
    if(numel(cl) ~= numel(cu)), error('cl and cu must be the same length'); end    
    %Add Objective
    try
        if(isa(con,'sym'))
            B.AddCon(con,cl,cu);
        else
            B.AddCon(con(x),cl,cu);
        end
    catch ME
        error('There was an error processing your constraints into Symbolic Form. Please examine the error below and correct your function:\n\nError: %s',ME.message);
    end
end

%Check Bounds & Integer Constraints
if(~isempty(lb) && numel(x0) ~= numel(lb)), error('x0 and lb must be the same length'); end
if(~isempty(ub) && numel(x0) ~= numel(ub)), error('x0 and ub must be the same length'); end
if(~isempty(xtype) && numel(x0) ~= numel(xtype)), error('x0 and xtype must be the same length'); end
if(isempty(lb)), lb = -Inf(ndec,1); end
if(isempty(ub)), ub = Inf(ndec,1); end
B.AddBounds(x,lb,ub);
if(~isempty(xtype)), B.AddIntegers(x,xtype); end

%Build & Return Object
B.Build();
Opt = GetOPTI(B,sopts);
SymB = B;