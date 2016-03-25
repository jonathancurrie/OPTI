%% IPOPT PARFOR testing
clc
clear all

%Objective
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
%Constraints
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];
%Setup Options
opts = optiset('solver','ipopt');
Opt = opti('fun',fun,'bounds',lb,ub,'opts',opts);

%Generate a series of starting points
n = 1000;
X0 = [linspace(lb(1),lb(end),n); linspace(ub(1),ub(end),n)];

%Preallocate solution vector
sols = zeros(2,n);
solp = zeros(2,n);

%% run serial (8.17s)
tic
for i = 1:n
    sols(:,i) = solve(Opt,X0(:,i));
end
toc

%%
parpool(4);

%% run parfor (2.6s)
tic
parfor i = 1:n
    solp(:,i) = solve(Opt,X0(:,i));
end
toc

%% check results
err = norm(sols-solp)