%% Test Problems for DSDP
clc
clear

%% DSDP Options Function
clc
dsdpset

a = dsdpset

%% LP1 [-31.4]
clc
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
b = [16;28;6];    
lb = [0;0];
ub = [10;10];

[x,p,d,e,i] = dsdp(-f,A,b,lb,ub)

%% LP2 [2]
clc
f = [8,1]';
A = sparse([-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]); 
b = [-4,2,21,39,3,0,0]';

[x,p,d,e,i] = dsdp(-f,A,b)

%% LP3 [-3.75]
clc
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
b = [5, 5]';

[x,p,d,e,i] = dsdp(-f,A,b)

%% LP4 [-97.5] (won't solve)
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

[x,fval,e,i,X] = dsdp(-f,A,b,lb,ub)

% [x,ff] = clp(f,[A;Aeq],[-Inf(size(b));beq],[b;beq],lb,ub)

%% SDP1 [1.4142]
clc
clear
%Objective
f = 1;
%SDP [x sqrt(2); sqrt(2) x] >= 0
ind = triu(ones(2))==1;
A = eye(2);
C = [0 sqrt(2); sqrt(2) 0];
sdp = sparse([C(ind) -A(ind)])
%Options
opts.display=2;
[x,ff,e,i,X] = dsdp(-f,[],[],[],[],sdp,[],opts)

%% SDP2 [4]
clc
clear
%Objective
f = [1;1];
%Linear Constraints
lb = [0;0];
ub = [10;10];
%SDP Constraints [x1 2; 2 x2] >= 0
ind = triu(ones(2))==1;
C = [0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
sdp = sparse([C(ind) -A0(ind) -A1(ind)]);
%Setup Options
opts.display=2;
[x,ff,e,i,X] = dsdp(-f,[],[],lb,ub,sdp,[],opts)

%% SDP3 [1.2] (Note assume matrices are symmetric triu to avoid x5)
clc
clear
%Objective
f = [1 0 0 0]';
%SDP Constraint1 [x2 x3;0 x4] <= x1*eye(2)
ind = triu(ones(2))==1;
C = zeros(2);
A0 = eye(2);
A1 = -[1 0; 0 0];
A2 = -[0 1; 0 0];
A3 = -[0 0; 0 1];
sdp{1} = sparse([C(ind) -A0(ind) -A1(ind) -A2(ind) -A3(ind)]);
%SDP Constraint2 [x2 x3; 0 x4] >= [1 0.2; 0.2 1]
C = -[1 0.2; 0.2 1];
A0 = zeros(2);
A1 = [1 0; 0 0];
A2 = [0 1; 0 0];
A3 = [0 0; 0 1];
sdp{2} = sparse([C(ind) -A0(ind) -A1(ind) -A2(ind) -A3(ind)]);
%Setup Options
opts.display=2;
opts.zbar = 0;
opts.ptol = 1e-8;
opts.rho = 3;
[x,ff,e,i,X] = dsdp(-f,[],[],[],[],sdp,[],opts)

%% SDP4 [10.1787]
clc
clear
%Objective
f = [1 1 1];
%Linear Constraints
lb = [10;0;0];
ub = [1000;1000;1000];
%SDP Constraints [x1 1 2; 1 x2 3; 2 3 100] >= 0
ind = triu(ones(3))==1;
C = [0 1 2; 1 0 3; 2 3 100];
A0 = [1 0 0; 0 0 0; 0 0 0];
A1 = [0 0 0; 0 1 0; 0 0 0];
A2 = zeros(3);
sdp = sparse([C(ind) -A0(ind) -A1(ind) -A2(ind)]);
%Setup Options
opts.display=2;
[x,ff,e,i,X] = dsdp(-f,[],[],lb,ub,sdp,[],opts)

%% SDP 5 [2], fixed var
clc
clear
%Objective
f = [1;1];
%Linear Constraints
lb = [0;0];
ub = [10;10];
%SDP Constraints [x1 2; 2 x2] >= 0
ind = triu(ones(2))==1;
C = [0 2; 2 0];
A0 = [1 0; 0 0];
A1 = [0 0; 0 1];
sdp = sparse([C(ind) -A0(ind) -A1(ind)]);
%Setup Options
opts.display=2;
opts.fixed = [1 3]; %variable 1, fixed to 3.0
[x,ff,e,i,X] = dsdp(-f,[],[],lb,ub,sdp,[],opts)


%% SDP 6 (testing memory leaks - use task manager to view matlab ram)
% clc
% clear all
% %Number of cones
% n = 30000;
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
% sdp = repmat({sparse([C(ind) -A0(ind) -A1(ind)])},n,1);
% %Setup Options
% opts.display=2;
% [x,ff,e,i,X] = dsdp(-f,[],[],lb,ub,sdp,[],opts);

