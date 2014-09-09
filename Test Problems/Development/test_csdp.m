%% Test Problems for DSDP
clc
clear

%% CSDP Options Function
clc
csdpset

a = csdpset

%% LP1 [-31.4]
clc
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
lb = [0;0];
ub = [10;10];

opts = [];
opts.display = 2;
opts.maxiter = 15;
% opts.writeprob = 'prob.dat-s';
% opts.writesol = 'sol.dat-s';
[x,ff,e,i,X] = csdp(f,A,b,lb,ub,[],[],opts)

%% LP2 [2]
clc
f = [8,1]';
A = sparse([-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]); 
b = [-4,2,21,39,3,0,0]';

[x,p,d,e] = csdp(f,A,b)

%% LP3 [-3.75]
clc
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
b = [5, 5]';

[x,p,d,e] = csdp(f,A,b)

%% LP4 [-97.5]
clc
clear all
f = -[1 2 3]';
A = sparse([-1,1,1; 1,-3,1]);
b = [20,30]';
Aeq = sparse([1 1 1]);
A = [A;Aeq;-Aeq];
beq = 40;
b = [b;beq;-beq];
lb = [0 0 0]';
ub =[40 inf inf]';

[x,fval,e,i] = csdp(f,A,b,lb,ub)

% [x,ff] = clp(f,[A;Aeq],[-Inf(size(b));beq],[b;beq],lb,ub)

%% SDP1 [1.4142]
clc
clear
%Objective
f = 1;
%SDP [x sqrt(2); sqrt(2) x] >= 0
A = eye(2);
C = -[0 sqrt(2); sqrt(2) 0];
sdp = sparse([C(:) A(:)]);
%Options
opts.display=2;
[x,ff,e,i] = csdp(f,[],[],[],[],sdp,[],opts)

%% SDP 1b Johan [-1]
clc
clear
%Objective
f = 1;
%SDP [1 x; x 1] >= 0
A = [0 1; 1 0];
C = -eye(2);
sdp = sparse([C(:) A(:)]);
%Options
opts.display = 2;
tic
[x,ff,e,i] = csdp(f,[],[],[],[],sdp,[],opts)
toc

%% SDP2 [4]
clc
clear
%Objective
f = [1;1];
%Linear Constraints
lb = [0;0];
ub = [10;10];
%SDP Constraints [x1 2; 2 x2] >= 0
C = -[0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
sdp = sparse([C(:) A0(:) A1(:)]);
%Setup Options
opts.display=2;
[x,ff,e,i] = csdp(f,[],[],lb,ub,sdp,[],opts)

%% SDP3 [1.2] (Note assume matrices are symmetric triu to avoid x5)
clc
clear
%Objective
f = [1 0 0 0]';
%SDP Constraint1 [x2 x3;x3 x4] <= x1*eye(2)
C = zeros(2);
A0 = eye(2);
A1 = -[1 0; 0 0];
A2 = -[0 1; 1 0];
A3 = -[0 0; 0 1];
sdp{1} = sparse([C(:) A0(:) A1(:) A2(:) A3(:)]);
%SDP Constraint2 [x2 x3; x3 x4] >= [1 0.2; 0.2 1]
C = [1 0.2; 0.2 1];
A0 = zeros(2);
A1 = [1 0; 0 0];
A2 = [0 1; 1 0];
A3 = [0 0; 0 1];
sdp{2} = sparse([C(:) A0(:) A1(:) A2(:) A3(:)]);
%Setup Options
opts.display=2;
[x,ff,e,i] = csdp(f,[],[],[],[],sdp,[],opts)

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
sdp = sparse([C(:) A0(:) A1(:) A2(:)]);
%Setup Options
opts.display=2;
[x,ff,e,i] = csdp(f,[],[],lb,ub,sdp,[],opts)


%% SDP CSDP Example
clc
clear
%Objective
f = -[1 2];
%SDP Constraint 1
C = [2 1; 1 2];
A0 = -[3 1; 1 3];
A1 = -zeros(2);
sdp{1} = sparse([C(:) A0(:) A1(:)]);
%SDP Constraint 2
C = [3 0 1; 0 2 0; 1 0 3];
A0 = -zeros(3);
A1 = -[3 0 1; 0 4 0; 1 0 5];
sdp{2} = sparse([C(:) A0(:) A1(:)]);
%SDP Constraint 3
C = zeros(2);
A0 = -[1 0; 0 0];
A1 = -[0 0; 0 1];
sdp{3} = sparse([C(:) A0(:) A1(:)]);
%Setup Options
opts.display=2;
% opts.writesol='sol.dat-s';
[x,ff,e,i,X] = csdp(f,[],[],[],[],sdp,[],opts)

X{1}
X{2}
X{3}

%% SDP5
% clc
% clear
% load dsdpdebug
% 
% opts.display=2;
% [y,fvals,exitflag,stats,X] = csdp(model.f,[],[],[],[],model.sdcone,[],opts)

%% SDP 6 (testing memory leaks - use task manager to view matlab ram)
% clc
% clear all
% %Number of cones
% n = 10000;
% %Objective
% f = [1;1];
% %Linear Constraints
% lb = [0;0];
% ub = [10;10];
% %SDP Constraints [x1 2; 2 x2] >= 0
% ind = triu(ones(2))==1;
% C = [0 2; 2 0];
% A0 = [1 0; 0 0];
% A1 = [0 0; 0 1];
% sdp = repmat({sparse([C(:) A0(:) A1(:)])},n,1);
% %Setup Options
% opts.display=2;
% opts.maxtime = 10;
% [x,ff,e,i,X] = csdp(f,[],[],lb,ub,sdp,[],opts);

%% YALMIP comparisons
% %%
% x = sdpvar(1,1);
% %y = sdpvar(1,1);
% z = sdpvar(1,1);
% 
% X = [x 2;2 z];
% 
% F = set('X>=0');
% % F = F+set('x>=1');
% F = F+set('z>=1');
% F = F+set('x<=10');
% F = F+set('z<=10');
% sol = solvesdp(F,x+z)
% 
% %%
% [At,b1,c1,K]=readsdpa('complete2uub.dat-s')
% 
% A1 = full(At)
% 
% %%
% [At,b2,c2,K]=readsdpa('prob.dat-s')
% 
% A2 = full(At)
% 
% 
% %%
% t = sdpvar(1,1);
% Y = sdpvar(2,2);
% F = set('Y<=t*eye(2)');
% F = F+set('Y>=[1 0.2;0.2 1]');
% sol = solvesdp(F,t)
% 
% %%
% [At,b1,c1,K]=readsdpa('test.dat-s')
% 
% A1 = full(At)
% 
% %%
% [At,b2,c2,K]=readsdpa('prob.dat-s')
% 
% A2 = full(At)
% 
% 
% %%
% clc
% x = sdpvar(1,1);
% y = sdpvar(1,1);
% 
% con = [x + 4*y <= 16; 6*x + 4*y <= 28, 2*x - 5*y <= 6; 0 <= x <= 10; 0 <= y <= 10];
% obj = -6*x - 5*y;
% 
% opts = sdpsettings('solver','csdp');
% sol = solvesdp(con,obj,opts)
% 
% double(x)
% double(y)
% double(obj)
% 
% %%
% [At,b1,c1,K]=readsdpa('lp1.dat-s')
% 
% A1 = full(At)
% 
% %%
% [At,b2,c2,K]=readsdpa('prob.dat-s')
% 
% A2 = full(At)
