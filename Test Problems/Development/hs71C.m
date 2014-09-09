function [c ceq gradc gradceq] = hs71C(x)

c = -(prod(x)-25);
ceq = sum(x.^2)-40;

if(nargout > 2)
    gradc = -(prod(x)./x')';
    gradceq = 2*x;
end
