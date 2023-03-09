function [prob,p] = amplRead(filename,extraArgs,amplPath,isNLP)
%AMPLREAD  Read an AMPL File and converts it to Matlab Values / Functions
%
% prob = amplRead(filename) reads the .nl or .mod file specified by filename 
% and converts it to an optiprob. If the file is a .mod the AMPL executable
% is used to convert it to a .nl, before reading it into optiprob. If the
% AMPL executable cannot be found on the MATLAB path, you will need to
% enter its path as the third argument (see below). The returned structure 
% is solver independent so the user can manually extract matrices as 
% required or supply it directly to opti().
%
% prob = amplRead(filename,extraArgs) allows you to pass extra arguments to 
% AMPL when converting your model file from .mod to .nl, such as options or
% data. Pass multiple arguments as a cell array. Files not on the MATLAB
% path must be passed as a full path.
%
% prob = amplRead(filename,extraArgs,amplPath) allows you to specify the
% full path to the AMPL executable on your system.
%
% prob = amplRead(filename,extraArgs,amplPath,isNLP) allows you to specify
% if the problem should be read as a (MI)NLP, even if identified by the
% asl interface as (MI)LP/QP/QCQP.
%
% The routine underneath uses the Netlib AMPL Solver Library code. You must 
% have a licensed version of AMPL (or the student edition) present on your 
% computer to read .mod files. For more information see www.ampl.com.

%   Copyright (C) 2011-2013 Jonathan Currie (Control Engineering)

%Optional args
if(nargin < 4), isNLP = 0; end
if(nargin < 3), amplPath = []; end
if(nargin < 2), extraArgs = []; end

if(~ischar(filename))
    error('Filename must be a char array!');
end
  
% Check if we are processing a .mod file
if(~isempty(strfind(filename,'.mod')))
    filename = amplConvert(filename,extraArgs,amplPath);
end

% Make sure we have file extension (always assume .nl)
if(isempty(strfind(filename,'.nl')))
    filename = [filename '.nl'];
end 

% If filename includes absolute path, we can skip which
if(~isempty(strfind(filename,':')))
    p = filename;
    %Check it exists
    if(~exist(filename,'file'))
        error('Cannot locate file %s!',filename);
    end
else %Locate the full path to the file (must be on matlab path)
    p = which(filename);
    %Check it exists
    if(isempty(p))
        error('Cannot locate file %s!',filename);
    end
end
if(isNLP), isNLP = 1; else isNLP = 0; end

%Read in Problem Data
[aslprob,sizes] = asl('open',p,isNLP);
%Assign args
x0 = aslprob.x0;
lb = aslprob.lb;
ub = aslprob.ub;
cl = aslprob.cl;
cu = aslprob.cu;
sense = aslprob.sense;
objbias = aslprob.objbias;

%Deallocate Sizes Vector
ndec    = sizes(1); %decision vars
ncon    = sizes(2); %constraints
% nzJac   = sizes(3); %non zeros in jacobian
lnc     = sizes(4); %linear network
nbv     = sizes(5); %linear binary vars
niv     = sizes(6); %linear integer vars
nlc     = sizes(7); %nonlinear constraints
nlnc    = sizes(8); %nonlinear network constraints
nlo     = sizes(9); %nonlinear objective
nlvb    = sizes(10); %nonlinear in both objective and constraints
nlvc    = sizes(11); %nonlinear vars in cons
nlvo    = sizes(12); %nonlinear vars in obj
nlvbi   = sizes(13); %integer nonlinear in both objective and constraints
nlvci   = sizes(14); %integer nonlinear just in constraints
nlvoi   = sizes(15); %integer nonlinear just in objectives
nwv     = sizes(16); %no lin arcs

if(nwv || lnc || nlnc)
    asl('close');
    error('Currently the OPTI AMPL Interface does not support network constraints or linear arcs');
end

%Process Integer Constraints
[~,~,int] = amplVars(ndec,nlvc,nlvo,niv,nbv,nwv,nlvb,nlvbi,nlvci,nlvoi);

% -- Construct OPTI Problem based on above sizes --%
doClose = 1;
%Check for LP
if(ncon && ~nlo && ~nlc && ~isNLP)
    %Build return problem
    prob = optiprob('f',aslprob.f,'lin',aslprob.A,cl,cu,'bounds',lb,ub,'x0',x0,'sense',sense,'objbias',objbias,'int',int);
    
%Check for QP
elseif(~isempty(aslprob.H) && isempty(aslprob.qcind) && ~isNLP)
    %Build return problem
    prob = optiprob('H',aslprob.H,'f',aslprob.f,'lin',aslprob.A,cl,cu,'bounds',lb,ub,'x0',x0,'sense',sense,'objbias',objbias,'int',int);  
    
%Check for QCQP
elseif(~isempty(aslprob.H) && ~isempty(aslprob.qcind) && ~isNLP)
    %Generate quadratic constraints
    [Q,l,qrl,qru,ind] = genQC(aslprob);
    %Remove quadratic constraints from linear constraints
    aslprob.A(ind,:) = [];
    cl(ind) = []; cu(ind) = [];
    %Build return problem
    prob = optiprob('H',aslprob.H,'f',aslprob.f,'lin',aslprob.A,cl,cu,'qcrow',Q,l,qrl,qru,'bounds',lb,ub,'x0',x0,'sense',sense,'objbias',objbias,'int',int); 
    
%Unconstrained NLP
elseif(~ncon && (isempty(lb) || all(isinf(lb))) && (isempty(ub) || all(isinf(ub)))) 
    %Assign callbacks
    fun     = @(x) asl('fun',x);
    grad    = @(x) asl('grad',x);    
    %Build return problem
    prob = optiprob('fun',fun,'grad',grad,'x0',x0,'sense',sense,'name','OPTI AMPL Problem');
    doClose = 0;

%NLP
else 
    %Assign callbacks
    fun     = @(x) asl('fun',x);
    grad    = @(x) asl('grad',x);
    hess    = @(x,sigma,lambda) asl('hess',x,sigma,lambda,1); %get sparse hess
    hstr    = @() asl('hesstr'); 
    %Ensure not just bounded
    if(~isempty(cl))  
        %Check for all linear constraints
        if(all(aslprob.conlin==0) && ~isNLP)
            A = asl('jac',zeros(size(aslprob.x0)),1); %get sparse jac;
            rl = cl; cl = [];
            ru = cu; cu = [];
            jac = []; jacstr = []; con = [];
        else
            con     = @(x) asl('con',x);
            jac     = @(x) asl('jac',x,1); %get sparse jac
            jacstr  = @()  asl('jacstr');  %get jac structure
            A = []; rl = []; ru = [];      %no linear constraints [not worth splitting off linear ones if they exist I suspect]
        end
    else
        A = []; rl = []; ru = []; jac = []; jacstr = []; con = []; cl = []; cu = [];
    end
    
    %Build return problem
    prob = optiprob('fun',fun,'grad',grad,'H',hess,'Hstr',hstr,'lin',A,rl,ru,'nl',con,cl,cu,'nljac',jac,'nljacstr',jacstr,...
                    'bounds',lb,ub,'x0',x0,'sense',sense,'int',int,'name','OPTI AMPL Problem');
    doClose = 0;
end

%Only close interface if we have all the data (i.e. linear or quadratic)
if(doClose)
    %asl('close');    %LEAVE OPEN so we can write the solution file
end
%Add full path to .nl file to structure
prob.path = p;
%Add constraint linearity 
prob.conlin = aslprob.conlin;

%Convert AMPL variables to opti types
function [ind_nl,ind_nlcon,int] = amplVars(ndec,nlvc,nlvo,niv,nbv,nwv,nlvb,nlvbi,nlvci,nlvoi)

%Ordering of AMPL Vars
% nonlinear max(nlvc, nlvo); see Table 4.
% linear arcs nwv
% other linear n_var - (max {nlvc, nlvo} + niv + nbv + nwv)
% binary nbv
% other integer niv

%Base integer string
int = repmat('C',1,ndec);

%Sort out #s of vars
nnl = max(nlvc,nlvo);
nla = nwv;
nln = ndec - (max (nlvc, nlvo) + niv + nbv + nwv);

%NL Indices
ind_nl  = 1:nnl;
%Allocate Linear Binary / Integer vars
i = nnl+nla+nln+1;
if(nbv)
    int(i:i+nbv-1) = 'B';
    i = i + nbv;
end
if(niv)
    int(i:i+niv-1) = 'I';
end

% Ordering of AMPL NL Vars
% continuous in an objective and in a constraint nlvb - nlvbi
% integer in an objective and in a constraint nlvbi
% continuous just in constraints nlvc - (nlvb + nlvci)
% integer just in constraints nlvci
% continuous just in objectives nlvo - (nlvc + nlvoi)
% integer just in objectives nlvoi

%Sort out #s of nl vars
nnlbthc = nlvb-nlvbi;
nnlbthi = nlvbi;
nnlcc = nlvc - (nlvb + nlvci);
nnlci = nlvci;
nnloc = nlvo - (nlvc + nlvoi);
nnloi = nlvoi;

%Determine nl constraints
ind_nlcon = 1:(nnlbthc + nnlbthi + nnlcc + nnlci);

%Allocate Nonlinear Integer Vars
i = nnlbthc + 1;
if(nnlbthi)
    int(i:i+nlvbi-1) = 'I';
    i = i + nnlbthi;
end
i = i + nnlcc;
if(nnlci)
    int(i:i+nnlci-1) = 'I';
    i = i + nnlci;
end
i = i + nnloc;
if(nnloi)
    int(i:i+nnloi-1) = 'I';
%     i = i + nnloi;
end


%Generate Quadratic Constraints from asl return structure
function [Q,l,qrl,qru,ind] = genQC(prob)
%Get Indices
ind = prob.qcind;
cl = prob.cl;
cu = prob.cu;
%Check for single constraint
if(length(ind) == 1)
    Q = prob.Q{1};
    l = prob.l;                
    qrl = cl(ind); 
    qru = cu(ind);
%Multiple constraints
else
    %Determine no of resulting QCs
    nqc = length(ind); 
    %Preallocate
    Q = cell(nqc,1); l = zeros(size(prob.l,1),nqc); 
    qrl = zeros(1,nqc); qru = zeros(1,nqc);
    for i = 1:length(ind)
        qci = ind(i);  
        Q(i) = prob.Q(i); 
        l(:,i) = prob.l(:,i);                
        qrl(i) = cl(qci);
        qru(i) = cu(qci);                      
    end
end







