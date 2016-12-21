function [f,g] = detGrad(fun,x0,xdata)
%DETGRAD  Determine whether supplied function contains gradient information

%   Copyright (C) 2011 Jonathan Currie (IPL)

%Sort out gradient (unfortunately rather inefficient)
no = nargout(fun);

if(no == 1)
    f = fun;
    g = [];
elseif(no == 2)
        f = @(x) fval(fun,x);
        g = @(x) gval(fun,x);
elseif(no < 0)
    %Either anonymous function or varargout (i.e. deal?)
    try
        if(nargin(fun) > 1)
            fun(x0,xdata);
        else
            fun(x0);       
        end
        f = fun;
        g = [];
    catch me
        try
            [~,~] = fun(x0);
            f = @(x) fval(fun,x);
            g = @(x) gval(fun,x);
        catch 
            error(['Unknown function / gradient format: ' me.message])
        end
    end
    
else
    error('Unknown function / gradient format');
end

%Dummy functions
function f = fval(fun,x)
[f,~] = fun(x);

function g = gval(fun,x)
[~,g] = fun(x);
