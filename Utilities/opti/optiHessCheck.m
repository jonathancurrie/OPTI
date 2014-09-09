function ok = optiHessCheck(Hess,grad,Jac,x0,name,warn,pattern)
%OPTIHESSCHECK  Check User Supplied Hessian of the Lagrangian and Optionally Pattern against mklJac
%
%   ok = optiHessCheck(Hess,grad,Jac,x0,name,warn,pattern)

%Tolerance
tol = 1e-3;

%Default input args and checks
if(nargin < 7), pattern = []; end
if(nargin < 6), warn = 1; end
if(nargin < 5), name = ''; end
if(nargin < 4), error('You must supply at least four arguments to this function'); end

%Default return
ok = true;

%Check user arguments
if(~isa(Hess,'function_handle')), derError(name,'This function requires function handle inputs'); end
%Hessian may have one or three inputs
switch(nargin(Hess))
    case 1
        Hess = @(x,sigma,lambda) Hess(x);
    case 3
        %all ok
    otherwise
        derError(name,'This function requires the user Hessian to one inputs (x) or three inputs (x,sigma,lambda)');
end
if(~isa(grad,'function_handle')), derError(name,'This function requires function handle inputs'); end
if(nargin(grad) ~= 1), derError(name,'This function requires the user gradient to have a single input argument (x)'); end
if(~isempty(Jac) && ~isa(Jac,'function_handle')), derError(name,'This function requires function handle inputs'); end
if(~isempty(Jac) && nargin(Jac) ~= 1), derError(name,'This function requires the user Jacobian to have a single input argument (x)'); end
if(isempty(x0) || issparse(x0) || ~isreal(x0) || ~isa(x0,'double') || (size(x0,1) > 1 && size(x0,2) > 1))
    derError('x0 must be a real, dense, double precision vector');
end
if(~isempty(pattern) && (~issparse(pattern) || ~isnumeric(pattern)))
    derError(name,'The sparsity pattern to be a sparse matrix');
end

%Default Warning String
if(isempty(name))
    wstr = sprintf('OPTI Derivative Checker detected a problem in a user supplied derivative\n');
else
    wstr = sprintf('OPTI Derivative Checker detected a problem in ''%s''\n',name);
end
isOK = true; %just used for warnings

%If we have nljac, then we can determine length of required lambda
if(~isempty(Jac))
    ju = Jac(x0);
    %Number of rows should = ncon
    lambda = ones(size(ju,1),1);
    lambda0 = zeros(size(ju,1),1);
else
    lambda = [];
    lambda0 = [];
end

%Evaluate full Hessian, check for size, NaN, Inf, keep for use later in
%pattern check
Huser = Hess(x0,1,lambda);
if(any(any(isnan(Huser))) || any(any(isinf(Huser))))
    derError(name,sprintf('One or more elements in the user supplied Hessian of the Lagrangian are NaN or Inf.\n\nPlease correct your initial guess to ensure a valid estimation of derivatives.'));
end
[rhu,chu] = size(Huser);
if(rhu~=chu), derError(name,'The user supplied Hessian of the Lagrangian is not square!'); end

%Decide we have a tril Hessian
isTril = false;
if(any(any(triu(Huser,1)))) %elements in upper tri
    isTril = true;
end

%Evaluate User Hessian just at objective
Huser_obj = Hess(x0,1,lambda0);

%Evaluate 2nd Derivative of Objective
if(~isempty(strfind(char(grad),'mklJac')))
    if(warn), optiwarn('OPTI:DoubleFiniteDiff','The derivative checker does not function correctly when the objective gradient uses mklJac - Results may be inaccurate.');end
    Hest_obj = mklJac(grad,x0,[],1e-3); %low precision if double finite difference
else
    Hest_obj = mklJac(grad,x0);
end

%Check Size of Estimated Objective Hessian
[rhobj,chobj] = size(Hest_obj);
if(rhobj~=chobj), derError(name,'The estimated Hessian of the Lagrangian is not square!'); end
if(rhobj~=rhu), derError(name,'The estimated Hessian of the Lagrangian does not have the same number of rows as the user supplied Hessian'); end
    
%Default warning String
wstr = sprintf('%s - There appear to be errors in your Hessian, please check the following elements:\n',wstr);

%Perform error check of objective Hessian only
[isOK,ok,wstr] = checkHess(Huser_obj,Hest_obj,'Objective',0,isTril,tol,isOK,ok,wstr);


%If we Jac, check each derivative
if(~isempty(Jac))
    if(~isempty(strfind(char(Jac),'mklJac')))
        if(warn), optiwarn('OPTI:DoubleFiniteDiff','The derivative checker does not function correctly when the constraint Jacobian uses mklJac - Results may be inaccurate.');end
        prec = 1e-3; %low precision if double finite difference
    else
        prec = [];
    end    
    for i = 1:length(lambda)
        lam = lambda0;
        lam(i) = 1;
        %Evaluate User Supplied
        Huser_con = Hess(x0,0,lam);
        %Evaluate Estimate
        Hest_con = hessJacEval(Jac,x0,i,prec); %slack - might be able to do the finite difference once
        %Check
        [isOK,ok,wstr] = checkHess(Huser_con,Hest_con,'Constraint',i,isTril,tol,isOK,ok,wstr);
    end                
end


%Check sparsity pattern if supplied
if(~isempty(pattern))
    %Check Size
    [rp,cp] = size(pattern);
    if(rp ~= rhu || cp ~= chu)
        derError('The size of the Hessian Sparsity Structure does not match the size of the returned Jacobian!');
    end
    Hstr = double(Huser ~= 0); 
    %If we subtract the estimated pattern from the original pattern we
    %should only have 0s (element found / never existed), and 1s (element
    %not identified - to be expected). But if we have any -1s, then we found
    %a new element not in the original pattern, error!
    newElems = (double(pattern) - Hstr) == -1;
    if(any(any((newElems))))
        isOK = false;
        %Build warning string
        wstr = sprintf('%s - There appear to be errors in your Hessian Sparsity Structure, OPTI identified elements in the following locations:\n',wstr);
        [I,J] = find(newElems);
        for i = 1:sum(sum(newElems))
            wstr = sprintf('%s    HessStruc(%d,%d) should equal 1\n',wstr,I(i),J(i));
        end
        ok = false;
    end
end

%Display Warning
if(~isOK && warn), optiwarn('OPTI:IncorrectHess',wstr); end
if(isOK), optiinfo('OPTI Derivative Checker detected no problems in ''%s''',name); end
    
function derError(name,msg)
throwAsCaller(MException('OPTI:DERCHECK','OPTI Derivative Checker detected a problem in ''%s'':\n\n%s',name,msg));


function [isOK,ok,wstr] = checkHess(user,estim,name,no,isTril,tol,isOK,ok,wstr)
%Check user hessian vs estimated hessian

if(any(any(isnan(user))) || any(any(isinf(user))))
    derError(name,sprintf('One or more elements in the user Hessian of the Lagrangian (%s %d) are NaN or Inf.\n\nPlease correct your initial guess to ensure a valid estimation of derivatives.',name,no));
end
if(any(any(isnan(estim))) || any(any(isinf(estim))))
    derError(name,sprintf('One or more elements in the estimated Hessian of the Lagrangian (%s %d) are NaN or Inf.\n\nPlease correct your initial guess to ensure a valid estimation of derivatives.',name,no));
end

%Guess if user has supplied only tril part of Hessian (should always be symmetric)
if(isTril) %elements in upper tri
    err = abs(estim-user);
else
    err = abs(tril(estim)-user); %our Hessian will be full + symmetric
end
%Indices of bad elements
ind = err > tol;
if(any(any(ind)))
    isOK = false;
    %Build warning string    
    [I,J] = find(ind);
    if(strcmp(name(1:3),'Obj'))
        sname = 'Hess[Obj]';
    else
        sname = sprintf('Hess[Con:%d]',no);
    end
    for i = 1:sum(sum(ind))
        wstr = sprintf('%s    %s(%d,%d) Error: %g\n',wstr,sname,I(i),J(i),err(I(i),J(i)));
    end
    ok = false;
end


function H = hessJacEval(jac,x,i,prec)
%Evaluate Hessian of Particular Jacobian
H = mklJac(@(x) full(jacEval(x,jac,i)),x,[],prec);

function J = jacEval(x,jac,i)
%Evaluate Jacobian and return just one row
Jall = jac(x);
J = Jall(i,:);

