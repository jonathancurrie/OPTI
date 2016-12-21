function O = GetOPTI(B,opts)
%Build OPTI Object from SymBuilder Object
%
%   Called By SYMBUILDER Class

%   Copyright (C) 2012-2013 Jonathan Currie (IPL)

if(nargin < 2), opts = symbset; else opts = symbset(opts); end

%Check object is built
if(~IsBuilt(B))
    error('Please build the object using Draft() or Build()');
end

%Check we have an objective (no good without!)
if(B.noObjs == 0)
    error('You can only construct an OPTI object when you have specified an objective!');
elseif(B.noObjs > 1)
    error('You can only construct an OPTI object with a single objective');
end

%Clear existing callback function from memory
clear symb_cb symb_ccb

%Separate objective and constraint indices
olin = B.objlin(B.indobj);
clin = B.conlin; clin(B.indobj) = [];

%Default options
if(~isempty(opts.optiOpts))
    Oopts = optiset(opts.optiOpts,'warnings','off'); %overrides solver and display settings
    %check if user has specified a solver/display
    if(~strcmpi(opts.solver,'auto'))
        optiwarn('OPTI:SymBSolver','As you have supplied optiOpts to symbset, your solver selection within symbset has been overwritten by the solver specified in optiOpts.');
    end
    if(~strcmpi(opts.display,'off'))
        optiwarn('OPTI:SymBDisplay','As you have supplied optiOpts to symbset, your display setting within symbset has been overwritten by the display setting specified in optiOpts.');
    end
else
    Oopts = optiset(opts.optiOpts,'warnings','off','solver',opts.solver,'display',opts.display);
end

%Check for Linear Problem
if(olin == 1 && all(clin == 1))
    %Build OPTI object
    O = opti(GetLinProb(B),'options',Oopts);
%Check for Quadratic Pprogram    
elseif(olin == 2 && all(clin == 1))
    %Build OPTI object
    O = opti(GetQuadProb(B),'options',Oopts);
%Check for QCQP/QCLP
elseif(olin <= 2 && all(clin <= 2))   
    %Build OPTI object
    O = opti(GetQuadProb(B),'options',Oopts);
%Check for NLP
elseif(olin <= 3 && all(clin <= 3))
    %If we are building an MINLP, add BONMIN linear information
    if(any(B.xtype ~= 'C'))
        %Constraint Linearity
        cons_lin = clin;
        cons_lin(clin > 1) = 0;
        %Variable linearity
        varlin = ones(size(B.jac,2),1);
        nlvars = symvar(B.jac); %Variables in the jacobian are nonlinear WHAT ABOUT GRADIENT?
        for i = 1:length(nlvars)
            varlin(logical(nlvars(i) == B.vars)) = 0;
        end
        if(~isempty(Oopts.solverOpts))
            Oopts = optiset(Oopts,'warnings','off','solverOpts',bonminset(Oopts.solverOpts,'cons_lin',cons_lin,'var_lin',varlin));
        else
            Oopts = optiset(Oopts,'warnings','off','solverOpts',bonminset('cons_lin',cons_lin,'var_lin',varlin));
        end
    end
    %Build OPTI object
    O = opti(GetNLProb(B,opts),'options',Oopts);
else
    error('problem type not yet supported');
end