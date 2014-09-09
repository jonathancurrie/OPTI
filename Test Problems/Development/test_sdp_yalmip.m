%% SDP Test YALMIP
clc
clear

%% SDP1
clc
clear
x = sdpvar(1,1);
%SDP Constraints [x sqrt(2); sqrt(2) x] >= 0
X = [x sqrt(2); sqrt(2) x];
% Solve
opts = sdpsettings('solver','csdp','verbose',1);
sol = solvesdp(X>=0,x,opts)

double(X)


%% SDP2
clc
x = sdpvar(1,1);
z = sdpvar(1,1);

X = [x 2;2 z];
F = set('X>=0');
F = F+set('x>=0');
F = F+set('z>=0');
F = F+set('x<=10');
F = F+set('z<=10');
sol = solvesdp(F,x+z)

plot(F,[x;z])

%%
clc
x = sdpvar(1,1);
solvesdp([1 x(1);x(1) 1]>=0,sum(x(1)),sdpsettings('solver','csdp'))

%%






