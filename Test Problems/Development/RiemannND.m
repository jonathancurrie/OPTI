function R = RiemannND(x)
% Riemann's non-differentiable function, R(x) for rational x = p/q

[p,q] = rat(x);
R = 0;
for k = 1:q-1
    R = R + sin(k^2*p*pi/q)/(sin(k*pi/2/q))^2;
end
R = R*pi/4/q/q;
end