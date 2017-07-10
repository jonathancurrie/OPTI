%% 3D Plot Testing

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% LP2 [x0, unsolved]
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','x0',[7.3;1.16;26.163],'options',opts)
plot(Opt)

%% LP2 [x0, unsolved, multi scale]
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','x0',[7.3;1.16;26.163],'options',opts)
plot(Opt,[5 13 0 4 24 31])

%% LP2 [NO x0, unsolved]
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','options',opts)
plot(Opt)

%% LP2 [NO x0, unsolved, multi scale]
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','options',opts)
plot(Opt,[5 13 0 4 24 31])

%% LP2 [NO x0, unsolved, multi solve + multiplot]
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
%Setup Options
opts = optiset('solver','clp','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',[1 1.1 1.2],40,'bounds',[6 1 25]',[12 3 30]','options',opts)
multisolve(Opt)
multiplot(Opt)

%% QCQP1 [-0.4004]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [5;4];
Q = eye(3);
l = [2;2;2];
r = 0.5;
lb = [0;0;0];
ub = [40;inf;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('display','iter','solver','ipopt'))
[x,fval,exitflag,info] = solve(Opt)

plot(Opt)

%% MIQP2 [-2.75]
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','auto');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'int','CIC','options',opts)
try
    [x,fval,exitflag,info] = solve(Opt)
    plot(Opt)
catch
end

%% MIQCQP2 [-0.0942]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [5;4];
Q = eye(3);
l = [0;0;0];
r = 1;
lb = [0;0;0];
ub = [40;inf;inf];
int = 'CCI';
%Setup Options
opts = optiset('solver','auto','display','iter');
% Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'int',int,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

plot(Opt)


%% SDP4 [10.1787]
clc
clear
%Objective
f = [1 1 1];
%Linear Constraints
lb = [10;0;0];
ub = [1000;1000;1000];
%SDP Constraints [x1 1 2; 1 x2 3; 2 3 100] >= 0
C = -[0 1 2; 1 0 3; 2 3 100];
A0 = [1 0 0; 0 0 0; 0 0 0];
A1 = [0 0 0; 0 1 0; 0 0 0];
A2 = zeros(3);
sdcone = [C A0 A1 A2];
%Setup Options
dopts = dsdpset('ptol',1e-8,'rho',3,'zbar',0);
opts = optiset('solver','dsdp','display','iter','solverOpts',dopts);
% Solve
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% NLP1 Hock & Schittkowski #71
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','all','display','iter','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)
plot(Opt,2)

%% NLP1 Hock & Schittkowski #71 [multiplot]
clc
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','all','display','iter','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= multisolve(Opt,x0)
multiplot(Opt)
