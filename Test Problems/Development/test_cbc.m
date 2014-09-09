%% Test Problems for CBC
clc
clear

%% CBC Options Function
clc
cbcset

a = cbcset

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
[x,ff,e,i,c] = cbc([],f,A,rl,ru,lb,ub,[],[],[],opts)

%% LP Bounded
clc
f = -[6 5]';    
lb = [0;0];
ub = [10;10];

opts = [];
opts.display = 1;
opts.maxiter = 15;
[x,ff,e,i,c] = cbc([],f,[],[],[],lb,ub,[],[],[],opts)

%% LP Unbounded
clc
f = -[6 5]';   
A = sparse([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;6];  

opts = [];
opts.display = 1;
opts.maxiter = 15;
[x,ff,e,i,c] = cbc([],f,A,rl,ru,[],[],[],[],[],opts)


%% ILP1 [-29]
clc
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;6];    
lb = [0;0];
ub = [10;10];
xtype = 'II';

opts = [];
opts.display = 1;
opts.maxiter = 15;
[x,ff,e,i,c] = cbc([],f,A,rl,ru,lb,ub,xtype,[],[],opts);

%% ILP1 [-29] x0
clc
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;6];    
lb = [0;0];
ub = [10;10];
xtype = 'II';
x0 = [4;1];

opts = [];
opts.display = 2;
opts.maxiter = 15;
[x,ff,e,i,lambda] = cbc([],f,A,rl,ru,lb,ub,xtype,[],x0,opts);

%% LARGE ILP
clc
clear cbc
prob = coinRead('testMILP2.mps');

opts = cbcset('reduceAndSplitCuts','on','twoMirCuts','on');
opts.display = 1;
opts.maxtime = 30;
% opts.gomoryCuts = 'on';
% opts.maxnodes = 1e5;

[~,ff,e,i,lambda] = cbc(tril(prob.H),prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub,prob.int,[],[],opts)

%% LARGE ILP [no display crash testing]
clc
clear cbc
prob = coinRead('testMILP2.mps');

opts = [];%cbcset('reduceAndSplitCuts','on','twoMirCuts','on');
opts.display = 0;
opts.heuristics = 'On';
opts.gomoryCuts = 'On';
opts.preprocess = 'SOS';
opts.presolve = 'On';
opts.scaling = 'automatic';
opts.nodeStrategy = 'Fewest';
opts.strategy = 1;
opts.maxtime = 30;

[~,ff,e,i,cc] = cbc(tril(prob.H),prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub,prob.int,[],[],cbcset)

%% MILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];  
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% MILP2
clc
%Objective & Constraints
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0]; 
b = [20;30];  
Aeq = [0 1 0 -3.5];
beq = 0;
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'eq',Aeq,beq,'bounds',[0 0 0 2]',[40 inf inf 3]','int','CCCI','options',opts)
[x,fval,exitflag,info] = solve(Opt)
% plot(Opt)

%% MILP3 Infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = [4,1; 1.5,1; -2, 1; -0.2, 1]; 
b = [5.3;4.5;-2.5; 1.5];  
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% MILP4
clc
%Objective & Constraints
f = -[-1, 2]';
A = [2, 1;-4, 4];
b = [5, 5]';
e = -[1, 1];
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'mix',A,b,e,'int','II','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt)

%% MILP5
clc
%Objective & Constraints
f = [3, -7, -12]';
A = [-3, 6, 8;6, -3, 7;-6, 3, 3];
b = [12, 8, 5]';
e = [-1, -1, -1];
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'mix',A,b,e,'int','III','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% MILP6 
clc
%Objectie + Constraints
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];
%Special Ordered Sets
sos = '12';
sosind = {(1:2) (3:5)};
soswt = {(1:2) (3:5)};
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('f',f,'ineq',A,b,'ub',ub,'sos',sos,sosind,soswt,'options',opts)
[x,fval,exitflag,info] = solve(Opt)

% write(Opt,'testsos.mps')

%% BILP1
clc
%Objective & Constraints
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
b = [6;9;1;3];  
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BB','options',opts)
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt,3)

%% BILP2
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
b = [9; 1; 0; 0];
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('grad',f,'ineq',A,b,'int','BBBB','options',opts)
[x,fval,exitflag,info] = solve(Opt)

%% MIQP1 [-11.5]
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
b = [3; 3.5];
ivars = 1;
%Setup Options
opts = optiset('display','iter','solver','cbc');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'bounds',[0;0],[10;10],'int',ivars,'options',opts)
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,4)

%% MIQP2 [-2.75]
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
%Setup Options
opts = optiset('solver','cbc','display','iter');
%Build & Solve
Opt = opti('hess',H,'grad',f,'ineq',A,b,'int','CIC','options',opts)
[x,fval,exitflag,info] = solve(Opt)



%% Test MILP problems
clc
prob = milp_prob(10);
%Setup Options
copts = cbcset('useCuts','Off');
opts = optiset('solver','cbc','display','iter','solverOpts',copts);
%Build & Solve
Opt = opti(prob,opts)
[x,fval,exitflag,info] = solve(Opt)
