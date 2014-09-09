function [f,g] = hs71F(x)

f = x(1)*x(4)*sum(x(1:3)) + x(3);

if(nargout > 1)
    g = [ x(1)*x(4) + x(4)*sum(x(1:3))
                    x(1)*x(4)
                    x(1)*x(4) + 1
                    x(1)*sum(x(1:3)) ];
end