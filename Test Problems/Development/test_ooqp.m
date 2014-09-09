%% Testing OOQP
clc
clear

%% QP1 -2.83333333301227;
H = speye(3);
f = -[2 3 1]';
A = sparse([1 1 1;3 -2 -3; 1 -3 2]); 
b = [1;1;1];      

opts = [];
opts.display = 2;

[x,fval,e,i,l] = ooqp(H,f,A,-Inf(size(b)),b,[],[],[],[],opts)      

%% QP2 -8.22222220552525 
clc
H = (sparse([1 -1; -1 2]));
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1]);
b = [2; 2; 3]; 
lb = [0;0];

opts = [];
opts.display = 2;

[~,fval,e,i] = ooqp(tril(H),f,A,-Inf(size(b)),b,[],[],lb,[],opts)

%% QP3 -6.41379310344827
H = tril(sparse([1 -1; -1 2]));
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1]);
b = [2; 2; 3];
Aeq = sparse([1 1.5]);
beq = 2;
lb = [0;0];
ub = [10;10];  

opts = [];
opts.display = 2;
opts.linear_solver = 'pardiso';
[x,f,e,i] = ooqp(H,f,A,-Inf(size(b)),b,Aeq,beq,lb,ub,opts)

%% LP1 -3.75
clc
H = [];
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
b = [5, 5]';
   
opts = [];
opts.display = 2;
[x,f,e,i] = ooqp(H,f,A,-Inf(size(b)),b,[],[],[],[],opts) 

%% LP2 3
clc
H = [];
f = [1, 2, 3, 7, 8, 8];
A = -sparse([5, -3, 2, -3, -1, 2; -1, 0, 2, 1, 3, -3;1, 2, -1, 0, 5, -1]);
b = -[-5, -1, 3]';
lb = zeros(6,1);
ub = 10*ones(6,1);   

opts = [];
opts.display = 2;
[x,f,e,i,l] = ooqp(H,f,-A,-Inf(size(b)),-b,[],[],lb,ub,opts) 

%% LP3 4
clc
n = 40;
t = (0:n-1)';
y = 3.5 -.2*t;
b = y + 0.5*ones(size(y));
m = [ones(n,1),t(:)];
A = sparse([m,-m,eye(n)]);
H = [];spalloc(44,44,0);
f = [sum(m),sum(-m),2*ones(1,n)];
lb = zeros(n+4,1);
ub = [10, 10, 10, 10, 5*ones(1,n)];  

opts = [];
opts.display = 2;
[~,f,e,i,l] = ooqp(H,f,-A,-Inf(size(b)),-b,[],[],lb,ub,opts) 

%% LARGE LP
clc
prob = coinRead('maros-r7.mps');
opts = optiset('solver','ooqp','display','iter');
Opt = opti(prob,opts);
%Solve
[~,f,e,i] = solve(Opt)

%%
qp1 = qp_prob(1);
opts = optiset('solver','ooqp','display','iter');
Opt = opti(qp1,opts);
%Solve
[~,f,e,i] = solve(Opt)

