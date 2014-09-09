%% Test Problems for CLP
clc
clear

%% CLP Options Function
clc
clpset

a = clpset

%% LP1 [-31.4]
clc
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;6];    
lb = [0;0];
ub = [10;10];

opts = [];
opts.display = 1;
opts.maxiter = 15;
opts.algorithm = 'automatic';
opts.numThreads = 4;
[x,ff,e,i,lambda] = clp([],f,A,rl,ru,lb,ub,opts);

%% LP1 with redundant row [-31.4]
clc
f = -[6 5]';
A = sparse([1,4; 2 8; 6,4; 2, -5]); 
rl = -Inf(4,1);
ru = [16;32;28;6];    
lb = [0;0];
ub = [10;10];

opts = [];
opts.display = 1;
opts.maxiter = 15;
opts.algorithm = 3;
% opts.writeprob = 'prob.dat-s';
% opts.writesol = 'sol.dat-s';
[x,ff,e,i,lambda] = clp([],f,A,rl,ru,lb,ub,opts)

%% LP2 [2]
clc
f = [8,1]';
A = sparse([-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]); 
rl = -Inf(7,1);
ru = [-4,2,21,39,3,0,0]';

[x,ff,e,i,lambda] = clp([],f,A,rl,ru)

%% LP3 [-3.75]
clc
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
rl = -Inf(2,1);
ru = [5, 5]';

[x,ff,e,i,lambda] = clp([],f,A,rl,ru)

%% LP4 [-97.5]
clc
f = -[1 2 3]';
A = sparse([-1,1,1; 1,-3,1]);
b = [20,30]';
Aeq = sparse([1 1 1]);
A = [A;Aeq;-Aeq];
beq = 40;
ru = [b;beq;-beq];
rl = -Inf(size(ru));
lb = [0 0 0]';
ub =[40 inf inf]';

[x,ff,e,i,lambda] = clp([],f,A,rl,ru,lb,ub,opts)

%% LP5 [-30, no linear constraints]
clc
f = -[1, 2]';
lb = [0;0];
ub = [10;10];

[x,ff,e,i,lambda] = clp([],f,[],[],[],lb,ub,opts)


%% LP6 infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;-30];    
lb=[0;0];
ub=[10;10];

opts = [];
opts.display = 1;
opts.maxiter = 150;
opts.algorithm =5;

[x,ff,e,i,lambda] = clp([],f,A,rl,ru,lb,ub,opts)

%% LP7 unbounded
clc
f = -[1, 2]';
lb = [0;0];
ub = [10;Inf];

opts = [];
opts.display = 1;
opts.maxiter = 150;
opts.algorithm = 0;

[x,ff,e,i,lambda] = clp([],f,[],[],[],lb,ub,opts)

%% LARGE LP
clc
clear clp
prob = coinRead('maros-r7.mps');

opts = [];
opts.display = 1;
opts.maxiter = 10000;
opts.doPresolve = 1;
opts.algorithm = 5;
opts.numThreads = 4;

[~,ff,e,i,lambda] = clp(prob.H,prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub,opts)

%% AMPL Problem LP
clc
clear clp
prob = coinRead('prod.mps');

opts = [];
opts.display = 1;
opts.maxiter = 10000;
opts.algorithm = 0;
opts.doPresolve = 1;
opts.numThreads = 1;

[~,ff,e,i,lambda] = clp(prob.H,prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub,opts)

%% LP No Presolve Dual 
clc
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
rl = -Inf(2,1);
ru = [5, 5]';

opts = [];
opts.display = 1;
opts.maxiter = 15;
opts.algorithm = 0;
opts.doPresolve = 0;
opts.numThreads = 1;
[x,ff,e,i,lambda] = clp([],f,A,rl,ru,[],[],opts);

%% LP Empty after Presolve Abc
clc
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
rl = -Inf(2,1);
ru = [5, 5]';

opts = [];
opts.display = 1;
opts.algorithm = 0;
opts.numThreads = 2;
[x,ff,e,i,lambda] = clp([],f,A,rl,ru,[],[],opts);

%% QP1 -2.83333333301227;
clc
H = speye(3);
f = -[2 3 1]';
A = sparse([1 1 1;3 -2 -3; 1 -3 2]); 
b = [1;1;1];      

opts = [];
opts.display = 2;
opts.algorithm = 5;
opts.numThreads = 1;

% prob=optiprob('qp',H,f,'ineq',A,b);
% coinWrite(prob,'testqp.mps')

[x,fval,e,i,l] = clp(H,f,A,-Inf(size(b)),b,[],[],opts)      

%% QP2 -8.22222220552525 (w bias 1.7778)
clc
H = (sparse([1 -1; -1 2]));
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1]);
b = [2; 2; 3]; 
lb = [0;0];

opts = [];
opts.display = 1;
opts.algorithm = 3;
opts.objbias = 10;

[x,fval,e,i,l] = clp(tril(H),f,A,-Inf(size(b)),b,[],[],opts)  

%% QP3 -6.41379310344827
H = tril(sparse([1 -1; -1 2]));
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1]);
b = [2; 2; 3];
Aeq = sparse([1 1.5]);
beq = 2;
A = [A;Aeq]; rl = [zeros(size(b));beq]; ru = [b;beq];
lb = [0;0];
ub = [10;10];  

opts = [];
opts.display = 2;
[x,fval,e,i,l] = clp(H,f,A,rl,ru,lb,ub,opts)  
