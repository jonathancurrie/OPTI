%% Test Difference Routines

% Comparing Different Strategies

% 1) DERIVEST Suite
% 2) MKL djacobi
% 3) adiff

%% Test 1
clc
fun = @(x) 2*x(1) + 3*x(2);
x0 = [1;1];

tic
% g1 = gradest(fun,x0)
toc

tic
[g2,stat] = mklJac(fun,x0)
toc

tic
g3 = autoJac(fun,x0)
toc

%%
clc
f = symJac(fun)


%% Test 2
clc
fun = @(x) 2.1*x(1) + 3*x(2)^2 + 4.2*sqrt(x(3));
x0 = [1;1;2];

tic
% g1 = gradest(fun,x0)
toc

tic
[g2,stat] = mklJac(fun,x0)
toc

tic
g3 = autoJac(fun,x0)
toc

%% Test 3
clc
nlcon = @(x) [8 - x(1)^2 - x(2)^2 - x(3)^2 - x(4)^2 - x(1) + x(2) - x(3) + x(4);
              10 - x(1)^2 - 2*x(2)^2 - x(3)^2 - 2*x(4)^2 + x(1) + x(4);
              5 - 2*x(1)^2 - x(2)^2 - x(3)^2 - 2*x(1) + x(2) + x(4)];
x0 = [1;1;1;1];

tic
% g1 = jacobianest(nlcon,x0)
toc

tic
[jac,stat] = mklJac(nlcon,x0,3)
toc

tic
g3 = autoJac(nlcon,x0)
toc