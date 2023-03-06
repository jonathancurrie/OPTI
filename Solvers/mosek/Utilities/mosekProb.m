function prob = mosekProb(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,xint,sdcone,x0,opts)
%MOSEKPROB  From general inputs check then build MOSEK problem
%
%   prob = mosekProb(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,xint,sdcone,x0,opts)

%   This function is based in parts on examples from the MOSEK Toolbox, 
%   Copyright (c) 1998-2011 MOSEK ApS, Denmark.

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

%Get warning status
warn = strcmpi(opts.warnings,'all');

%Get sizes
ndec = length(f);
[nlin,Ad] = size(A);

%Error checks
if(~isempty(H))
    [r1,c1] = size(H);
    if(r1 ~= ndec || c1 ~= ndec), error('H is the wrong size, expected %d x %d',ndec,ndec); end
    %Symmetry check is done in OPTI, won't repeat here
end 
if(isempty(f)), error('You must supply the f vector to a MOSEK problem!'); end
if(xor(isempty(A),isempty(rl)) || xor(isempty(A),isempty(ru))), error('You must supply A + rl + ru for linear constraints'); end
if(~isempty(A))
    if(Ad ~= ndec), error('A is the wrong size, expected %d x %d',nineq,ndec); end
    if(nlin ~= length(rl)), error('The sizes of A and rl do not correspond'); end
    if(nlin ~= length(ru)), error('The sizes of A and ru do not correspond'); end
end
if(~isempty(lb) && length(lb) ~= ndec), error('lb is the wrong length, expected %d x 1',ndec); end
if(~isempty(ub) && length(ub) ~= ndec), error('ub is the wrong length, expected %d x 1',ndec); end
if(~isempty(xint))
    if(length(xint) ~= ndec), error('xint is the wrong length, expected %d x 1',ndec); end
    if(~ischar(xint)), error('xint should be a char array with characters ''C'', ''I'' or ''B'''); end    
end
if(~isempty(x0) && length(x0) ~= ndec), error('x0 is the wrong length, expected %d x 1',ndec); end

%Quadratic Objective
if(~isempty(H))
    if(~issparse(H))
        if(warn)
            optiwarn('mosek:sparse','The H matrix should be sparse, correcting: [sparse(H)]');
        end
        [prob.qosubi,prob.qosubj,prob.qoval] = find(sparse(tril(H)));    
    else
        [prob.qosubi,prob.qosubj,prob.qoval] = find(tril(H));
    end
end

%Linear Objective (if not sd)
if(isempty(sdcone))
    prob.c = f;
else
    prob.c = zeros(size(f));
end

%Linear constraints
prob.a = A;
if(isempty(prob.a))
    prob.a = sparse(0,length(f));
elseif(~issparse(prob.a))
    if(warn)
        optiwarn('mosek:sparse','The A matrix should be sparse, correcting: [sparse(A)]');
    end
    prob.a = sparse(prob.a);
end
if(size(f,2) > 1), f = f'; end
if(size(rl,2) > 1), rl = rl'; end
if(size(ru,2) > 1), ru = ru'; end
prob.blc = rl;
prob.buc = ru;

%Quadratic constraints
if(~isempty(Q))
    prob = genQC(prob,Q,l,qrl,qru,warn);
end

%Semidefinite Constraints
if(~isempty(sdcone))
    prob = genSDCONE(prob,f,sdcone);
end

%Integer variables
if(~isempty(xint))
    xint = lower(xint); [r,c] = size(xint);
    if(r > c), xint = xint'; end
    cind = strfind(xint,'c');
    bind = strfind(xint,'b');
    iind = strfind(xint,'i');
    if(length([cind bind iind]) ~= ndec) %double check no other chars
        error('xint should be a char array with characters ''C'', ''I'' or ''B''');
    end
    %Only need to setup if integer vars exist
    if(~isempty(bind) || ~isempty(iind))
        prob.ints.sub = sort([bind iind]); %sorted indices of all binary and integer vars
        %If binary, need to setup binary bounds
        if(~isempty(bind))
            if(isempty(lb)), lb = -Inf(size(f)); end %create false bounds if required
            lb(bind) = 0; 
            if(isempty(ub)), ub = Inf(size(f)); end %create false bounds if required
            ub(bind) = 1;
        end
    end
end

%Bounds (note must be done after binary stuff above)
prob.blx = lb;
prob.bux = ub;

%Initial guess
if(~isempty(x0))
    if(isfield(prob,'ints'))
        prob.sol.int.xx = x0;
    else
        prob.sol.itr.xx = x0;
    end
end

%Determine problem type
if(isempty(H)) %linear or sdp
    if(isempty(xint))
        if(isempty(sdcone))
            prob.names.name = 'OPTI MOSEK LP Problem';
        else
            prob.names.name = 'OPTI MOSEK SDP Problem';
        end
    elseif(all(lower(xint) == 'b'))
        prob.names.name = 'OPTI MOSEK BILP Problem';
    else
        prob.names.name = 'OPTI MOSEK MILP Problem';
    end
else %quadratic
    if(isempty(Q))
        if(isempty(xint))
            prob.names.name = 'OPTI MOSEK QP Problem';
        else
            prob.names.name = 'OPTI MOSEK MIQP Problem';
        end
    else %quadratically constrained
        if(isempty(xint))
            prob.names.name = 'OPTI MOSEK QCQP Problem';
        else
            prob.names.name = 'OPTI MOSEK MIQCQP Problem';
        end
    end
end


function prob = genQC(prob,Q,l,qrl,qru,warn)
%Generate MOSEK quadratic constraints

%Create default fields
prob.qcsubi = [];
prob.qcsubj = [];
prob.qcsubk = [];
prob.qcval = [];

ndec = length(prob.c);

%Error checks
if(~isempty(Q))
    if(isempty(l) || isempty(qrl) || isempty(qru))
        error('You must supply Q, l, qrl and qru for quadratic constraints');
    end
end

%Process cell by cell
if(iscell(Q))
   if(iscell(l) || iscell(qrl) || iscell(qru))
       error('Only Quadratic Constraints Q may be a cell!');
   end
   %Check overall sizes
   [r0,c0] = size(Q);
   [r1,c1] = size(l);
   [r2,c2] = size(qrl);
   [r3,c3] = size(qru);
   if(r0 > c0)
       Q = Q';
       [~,c0] = size(Q);
   end    
   if(r1 ~= c1) %not square
       if(r1 == c0)
           l = l';
           c1 = size(l,2);
       end
   end
   if(r2 > c2)
       qrl = qrl';
       [~,c2] = size(qrl);
   end
   if(r3 > c3)
       qru = qru';
       [~,c3] = size(qru);
   end
   if(c0 ~= c1 || c0 ~= c2 || c0 ~= c3)
       error('Quadratic Constraints Q + l + qrl + qru are not the same length!');
   end
   %For each cell
   for i = 1:c2
        Qt = Q{i};
        lt = l(:,i);
        rlt = qrl(i);
        rut = qru(i);
        %Error Check
        checkQC(Qt,lt,rlt,rut,ndec,sprintf('in cell %d',i),sprintf('in column %d',i));
        %Check sparsity
        if(~issparse(Qt))
            if(warn)
                optiwarn('mosek:sparse','The Q matrix in cell %d should be sparse, correcting: [sparse(Q)]',i);
            end
            Qt = sparse(Qt);
        end
        %Build structure fields
        [qcsubi,qcsubj,qcval] = find(tril(Qt));
        qcval = 2*qcval; %MOSEK includes 1/2x'Qx
        %Concatenate to original vectors
        prob.qcsubi = [prob.qcsubi; qcsubi];
        prob.qcsubj = [prob.qcsubj; qcsubj];
        prob.qcval = [prob.qcval; qcval];
        %Insert linear part into a matrix
        prob = insertQClin(prob,lt,rlt,rut,length(qcval));     
   end
else  %Single Quadratic Constraint
    if(size(l,2) > size(l,1))
        l = l';
    end       
    %Error Check
    checkQC(Q,l,qrl,qru,ndec,'','');
    %Check sparsity
    if(~issparse(Q))
        if(warn)
            optiwarn('mosek:sparse','The Q matrix should be sparse, correcting: [sparse(Q)]');
        end
        Q = sparse(Q);
    end
    %Build structure fields
    [prob.qcsubi,prob.qcsubj,prob.qcval] = find(tril(Q));
    prob.qcval = 2*prob.qcval; %MOSEK includes 1/2x'Qx
    %Insert linear part into a matrix
    prob = insertQClin(prob,l,qrl,qru,length(prob.qcval)); 
end

function prob = genSDCONE(prob,f,sdcone)
%Generate MOSEK semidefinite constraints
if(~iscell(sdcone)), sdcone = {sdcone}; end
ndec = length(f);
%Create default fields
prob.bardim    = [];
prob.barc.subj = []; %which cone
prob.barc.subk = []; %row
prob.barc.subl = []; %col
prob.barc.val  = [];
prob.bara.subi = []; %which A matrix
prob.bara.subj = []; %which cone
prob.bara.subk = []; %row
prob.bara.subl = []; %col
prob.bara.val  = [];

%For each sdcone
for i = 1:length(sdcone)
    dim = chkSDDim(sdcone{i},ndec,i);
    %Extract C, enter
    prob.bardim = [prob.bardim dim];
    [ci,cj,cv] = find(tril(reshape(sdcone{i}(:,1),dim,dim)));
    prob.barc.subj = [prob.barc.subj repmat(i,1,length(ci))];
    prob.barc.subk = [prob.barc.subk ci'];
    prob.barc.subl = [prob.barc.subl cj'];
    prob.barc.val = [prob.barc.val -cv'];
    %For each A, enter
    for j = 1:ndec
        [ai,aj,av] = find(tril(reshape(sdcone{i}(:,j+1),dim,dim)));
        prob.bara.subi = [prob.bara.subi repmat(j,1,length(ai))];
        prob.bara.subj = [prob.bara.subj repmat(i,1,length(ai))];
        prob.bara.subk = [prob.bara.subk ai'];
        prob.bara.subl = [prob.bara.subl aj'];
        prob.bara.val = [prob.bara.val av'];
    end        
end
%Modify a
prob.a = [prob.a; sparse(ndec,ndec)];
%Set bounds
prob.blc = [prob.blc;f];
prob.buc = [prob.buc;f];


function checkQC(Q,l,qrl,qru,ndec,msgC,msgV)
%Complete QC checks
if((ndec ~= size(Q,1)) || (ndec ~= size(Q,2)))
    error('Constraint Q matrix %s is the wrong size! Expected %d x %d',msgC,ndec,ndec);
end
if(ndec ~= size(l,1))
   error('Constraint l vector %s is the wrong size! Expected %d x 1',msgV,ndec);
end
if(~isscalar(qrl))
   error('Constraint qrl %s is the wrong size! Expected 1 x 1 ',msgV);
end
if(~isscalar(qru))
   error('Constraint qru %s is the wrong size! Expected 1 x 1 ',msgV);
end
try
    e = eig(Q);
catch
    e = eigs(Q);
end
if(any(e < 0))
   error('Constraint Q matrix %s is not positive semidefinite!',msgC);
end
if(all(e == 0))
   error('Constraint Q matrix %s is indefinite!',msgC);
end
sym = abs(tril(Q,-1)-triu(Q,1)') > 1e-9; %tol may need to be adjusted
if(any(any(sym)))       
    error('Constraint Q matrix %s is not symmetric',msgC);
end

function dim = chkSDDim(cone,ndec,i)
if(size(cone,2) ~= ndec+1)
    error('Semidefinite Cone %d does not have the correct number of columns',i);
end
dim = sqrt(size(cone,1));
if(floor(dim) ~= dim)
    error('Semidefinite Cone %d does not have the correct number of rows to form a square matrix');
end

function prob = insertQClin(prob,l,qrl,qru,n)
%Insert linear components of qc into prob structure

%Generate k vector
k = size(prob.a,1) + 1;
prob.qcsubk = [prob.qcsubk; k*ones(n,1)];
%Augment to linear a and constraint bounds
prob.a = [prob.a; l'];
prob.blc = [prob.blc; qrl];
prob.buc = [prob.buc; qru];

