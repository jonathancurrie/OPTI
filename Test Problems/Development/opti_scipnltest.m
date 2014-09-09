function res = opti_scipnltest(nlcon,x0)
% Testing SCIP Nonlinear Interface & SCIPVAR Instructions

if(isempty(x0))
    error('You must supply x0 to this function');
end

x = scipvar(size(x0));
n = nlcon(x);

if(length(n) > 1)
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
nl.cl = zeros(size(n));
nl.cu = ones(size(n));
nl.x0 = x0;
nl.nlcon_val = nlcon(x0);
if(any(isnan(nl.nlcon_val)) || any(isinf(nl.nlcon_val)))
    error('Starting guess results in NaN or Inf');
end

%Not used
f = ones(size(x0));
A = []; rl = []; ru = [];
lb = zeros(size(x0)); ub = ones(size(x0));
xtype = [];
qc =  [];
sos = [];
%Default Settings
opts.display = 0;
opts.testmode = 1;

res = scip([],f,A,rl,ru,lb,ub,xtype,sos,qc,nl,opts);