function prob = nrow2mix(prob,warn,doJac)
%NROW2MIX  Convert Nonlinear Row bounds to Mixed Nonlinear Constraints
%   prob = nrow2mix(prob)
%
%   e: -1 for <=, 0 for =, and 1 for >=

%   Copyright (C) 2012 Jonathan Currie (I2C2)

if(nargin < 3), doJac = true; end
if(nargin < 2), warn = 1; end    

%Assign common vars
nlcon = prob.nlcon;
cl = prob.cl;
cu = prob.cu;
nljac = prob.nljac;
nljacstr = prob.nljacstr;

if(xor(isempty(cl),isempty(cu)))
    error('You must supply both cl and cu!');
end
if(length(cl) ~= length(cu))
    error('cl and cu are not the same length!');
end
%Transpose as Required
if(size(cl,2) > 1)
    cl = cl';
end
if(size(cu,2) > 1)
    cu = cu';
end

%Indices
eq = cl == cu; neq = ~eq;
ile = isfinite(cu) & neq;
ige = isfinite(cl) & neq;

%If we have any dual row bounds, have to modify callbacks
ndbl = sum(ile&ige);
if(ndbl)
    %Create RHS + Type
    prob.nlrhs = [cu(ile); cl(ige); cl(eq)];
    on = ones(size(cl));
    prob.nle = [-1*on(ile); on(ige); 0*on(eq)];
    %New total number of constraints
    no = sum(ile+ige+eq); 
    if(isfield(prob,'sizes'))
        prob.sizes.ncon = prob.sizes.ncon + (no - (prob.sizes.nnlineq+prob.sizes.nnleq));
        prob.sizes.nnlineq = sum(ile+ige);
        prob.sizes.nnleq = sum(eq);
    end
    %Build new logical vectors
    fle = false(no,1);
    fge = false(no,1);
    feq = false(no,1);
    %Fill in indices
    fle(1:sum(ile)) = true; s = sum(ile);
    fge(s+(1:sum(ige))) = true; s = s + sum(ige);
    feq(s+(1:sum(eq))) = true;
    %Assign new constraint callback
    prob.nlcon = @(x) rowNlcon(x,nlcon,no,fle,ile,fge,ige,feq,eq);
    if(doJac)
        %Assign new jacstr callback
        if(~isempty(nljacstr))
            prob.nljacstr = @() rowNljacstr(nljacstr,no,fle,ile,fge,ige,feq,eq);
            nZ = nnz(prob.nljacstr());
        else
            nZ = []; %don't know
        end
        %Assign new jac callback
        if(~isempty(nljac))
            %If we are controlling Jac, then just reassign with the new no
            if(strcmp(char(nljac),'@(x) mklJac(prob.nlcon,x,nnl)'))
                prob.nljac = @(x) mklJac(nlcon,x,no);
            else
                prob.nljac = @(x) rowNljac(x,nljac,no,nZ,fle,ile,fge,ige,feq,eq);
            end
        end
        %Check for Hessian (haven't sorted out lambda yet...)
        if(~isempty(prob.H))
            if(warn), optiwarn('opti:dblbndHess',['Currently this interface does not support modifying lambda (for Hessian callbacks) '...
                                        'when converting "cl <= nlcon(x) <= cu" to "nlcon(x) nle nlrhs" with double bounded constraints '...
                                        '(both cl and cu different on one row). Removing Hessian Callback (sorry!)']); 
            end
            prob.H = [];
            prob.Hstr = [];
        end
    end
%Otherwise we can just get away with creating new vectors (much easier!)
else
    rhs = zeros(length(cl),1);
    e = zeros(length(cl),1);
    rhs(ile) = cu(ile); e(ile) = -1;
    rhs(ige) = cl(ige); e(ige) = 1;
    rhs(eq) = cl(eq);
    prob.nlrhs = rhs;
    prob.nle = e;
end

%Remove unused fields
prob.cl = []; prob.cu = [];

function con = rowNlcon(x,nlcon,no,fle,ile,fge,ige,feq,eq)
% Handle to allow conversion of row based nonlinear constraints to mix with double bounds
%Evaluate original function
c = nlcon(x);
%Assign return vector based on passed indices
con = zeros(no,1);
con(fle) = c(ile);
con(fge) = c(ige);
con(feq) = c(eq);


function jac = rowNljac(x,nljac,no,nZ,fle,ile,fge,ige,feq,eq)
% Handle to allow conversion of row based nonlinear jacobian to mix with double bounds
%Evaluate original function
c = nljac(x);
%Assign return vector based on passed indices
if(issparse(c) && ~isempty(nZ))
    jac = spalloc(no,length(x),nZ);
elseif(~sparse(c)) %ok if dense, if not, no preallocation
    jac = zeros(no,length(x));
end
jac(fle,:) = c(ile,:);
jac(fge,:) = c(ige,:);
jac(feq,:) = c(eq,:);


function str = rowNljacstr(nljacstr,no,fle,ile,fge,ige,feq,eq)
% Handle to allow conversion of row based nonlinear jacobian structure to mix with double bounds
%Evaluate original function
c = nljacstr();
%Assign return vector based on passed indices
if(~issparse(c))
    str = zeros(no,length(x));
end
%Won't bother with preallocating as only one call
str(fle,:) = c(ile,:);
str(fge,:) = c(ige,:);
str(feq,:) = c(eq,:);
