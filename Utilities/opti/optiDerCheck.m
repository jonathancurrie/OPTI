function ok = optiDerCheck(fun,grad,x0,name,warn,pattern)
%OPTIDERCHECK  Check User Supplied Derivatives and Optionally Pattern against mklJac
%
%   ok = optiDerCheck(fun,grad,x0,name,warn,pattern)

%Tolerance
tol = 1e-3;

%Default input args and checks
if(nargin < 6), pattern = []; end
if(nargin < 5), warn = 1; end
if(nargin < 4), name = ''; end
if(nargin < 3), error('You must supply at least three arguments to this function'); end

%Default return
ok = true;

%Check user arguments
if(~isa(fun,'function_handle')), derError(name,'This function requires function handle inputs'); end
if(nargin(fun) ~= 1), derError(name,'This function requires the user function to have a single input argument (x)'); end
if(~isa(grad,'function_handle')), derError(name,'This function requires function handle inputs'); end
if(nargin(grad) ~= 1), derError(name,'This function requires the user gradient to have a single input argument (x)'); end
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
isTrans = false;
isOK = true; %just used for warnings

%Evaluate Gradients
go = mklJac(fun,x0);
gu = grad(x0);

if(any(any(isnan(go))) || any(any(isinf(go))) || any(any(isnan(gu))) || any(any(isinf(gu))))
    derError(name,sprintf('One or more elements in the gradient are NaN or Inf.\n\nPlease correct your initial guess to ensure a valid estimation of derivatives.'));
end

% Ensure grad returns a double
checkIsDouble(gu, name);

%Check the same size
[ro,co] = size(go);
[ru,cu] = size(gu);
%Check for transposed
if(ro ~= co && ro == cu)
    wstr = sprintf('%s - It appears your gradient/Jacobian is transposed\n',wstr);
    gu = gu';
    isTrans = true; isOK = false;
    [ru,cu] = size(gu);
end
if(ro ~= co && co == ru)
    wstr = sprintf('%s - It appears your gradient/Jacobian is transposed\n',wstr);
    gu = gu';
    isTrans = true; isOK = false;
    [ru,cu] = size(gu);
end
if(ro ~= ru), derError(name,sprintf('Row size mismatch! OPTI #rows: %d, User #rows: %d',ro,ru)); end
if(co ~= cu), derError(name,sprintf('Column size mismatch! OPTI #columns: %d, User #columns: %d',co,cu)); end

%Check Abs Error
err = abs(go-gu); ind = err > tol;
if(any(any(ind)))
    isOK = false;
    %Build warning string
    wstr = sprintf('%s - There appear to be errors in your gradient/Jacobian, please check the following elements:\n',wstr);
    [I,J] = find(ind);
    for i = 1:sum(sum(ind))
        if(isTrans)
            wstr = sprintf('%s    Jac(%d,%d) Error: %g\n',wstr,J(i),I(i),err(I(i),J(i)));
        else
            wstr = sprintf('%s    Jac(%d,%d) Error: %g\n',wstr,I(i),J(i),err(I(i),J(i)));
        end
    end
    ok = false;
end

%Check sparsity pattern if supplied
if(~isempty(pattern))
    %Check size
    if(isTrans), pattern = pattern'; end
    [rp,cp] = size(pattern);
    if(rp ~= ru || cp ~= cu)
        derError(name,'The size of the Jacobian Sparsity Structure does not match the size of the returned Jacobian!');
    end    
    Jstr = double(gu ~= 0); 
    %If we subtract the estimated pattern from the original pattern we
    %should only have 0s (element found / never existed), and 1s (element
    %not identified - to be expected). But if we have any -1s, then we found
    %a new element not in the original pattern, error!
    newElems = (double(pattern) - Jstr) == -1;
    if(any(any((newElems))))
        isOK = false;
        %Build warning string
        wstr = sprintf('%s - There appear to be errors in your Jacobian Sparsity Structure, OPTI identified elements in the following locations:\n',wstr);
        [I,J] = find(newElems);
        for i = 1:sum(sum(newElems))
            wstr = sprintf('%s    JacStruc(%d,%d) should equal 1\n',wstr,I(i),J(i));
        end
        ok = false;
    end
end

%Display Warning
if(~isOK && warn), optiwarn('OPTI:IncorrectJac',wstr); end
if(isOK && warn), optiinfo('OPTI Derivative Checker detected no problems in ''%s''',name); end
    
function derError(name,msg)
throwAsCaller(MException('OPTI:DERCHECK','OPTI Derivative Checker detected a problem in ''%s'':\n\n%s',name,msg));

function checkIsDouble(v, name)
if (~isa(v, 'double'))
    if (isempty(name))
        error('OPTI:NotDouble', 'The %s derivative function does not return a double precision value (%s instead)', name, class(v));
    else
        error('OPTI:NotDouble', 'A user supplied derivative function does not return a double precision value (%s instead)', class(v));
    end
end