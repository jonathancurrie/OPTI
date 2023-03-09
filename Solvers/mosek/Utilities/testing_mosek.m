%% LP1 (x = [2.4 3.4], fval = -31.4)
clc
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
rl = -Inf(3,1);
ru = [16;28;6];   
lb = [0;0];
ub = [10;10];

[x,fval,exitflag,info] = moseklp(f,A,rl,ru,lb,ub)

%% LP2 (x = [10 2.5 27.5], fval = -97.5)
clc
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1; 1 1 1];
rl = [-Inf;-Inf;40];
ru = [20,30,40]';
lb = [0;0;0]';
ub = [40;inf;inf];
opts = mosekset('display','iter');

[x,fval,exitflag,info] = moseklp(f,A,rl,ru,lb,ub,[],opts)

%% LP3 (x = [0 2], fval = 2)
clc
f = [8,1]';
A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
rl = -Inf(7,1);
ru = [-4,2,21,39,3,0,0]';

[x,fval,exitflag,info] = moseklp(f,A,rl,ru)

%% BILP1 (x = [0 1], fval = -5)
clc
f = -[6 5]';
A = [-3,5; 6,4; 3, -5; -6, -4]; 
rl = -Inf(4,1);
ru = [6;9;1;3];  

[x,fval,exitflag,info] = mosekbilp(f,A,rl,ru)

%% BILP2 (x = [1 1 0 0], fval = -14)
clc
f = -[9 5 6 4]';
A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
rl = -Inf(4,1);
ru = [9; 1; 0; 0];

[x,fval,exitflag,info] = mosekbilp(f,A,rl,ru)


%% MILP1 (x = [4, 1], fval = -29)
clc
f = -[6 5]';
A = [1,4; 6,4; 2, -5]; 
rl = -Inf(3,1);
ru = [16;28;6];  
lb = [0;0];
ub = [10;10];
xint = 'II';

[x,fval,exitflag,info] = mosekmilp(f,A,rl,ru,lb,ub,xint)

%% MILP1a  (x = [1;2], fval = -3) (CHECK THIS IS WRONG IN 7.0.0.96?)
clc
f = -[-1, 2]';
A = [2, 1;-4, 4];
rl = -Inf(2,1);
ru = [5, 5]';
xint = 'II';
opts = mosekset('display','iter');

[x,fval,exitflag,info] = mosekmilp(f,A,rl,ru,[],[],xint,[],opts)

%% MILP2 (x = [40 10.5 19.5 3], fval = -122.5)
clc
f = -[1 2 3 1]'; 
A = [-1 1 1 10; 1 -3 1 0; 0 1 0 -3.5]; 
rl = [-Inf(2,1);0];
ru = [20;30;0];  
lb = [0;0;0;2];
ub = [40;inf;inf;3];
xint = 'CCCI';

[x,fval,exitflag,info] = mosekmilp(f,A,rl,ru,lb,ub,xint)

%% MILP3 Infeasible
clc
f = -[6 5]';
A = [4,1; 1.5,1; -2, 1; -0.2, 1]; 
rl = -Inf(4,1);
ru = [5.3;4.5;-2.5; 1.5];  
lb = [0;0];
ub = [10;10];
xint = 'II';

[x,fval,exitflag,info] = mosekmilp(f,A,rl,ru,lb,ub,xint)

%% MILP4 SOS (No SOS in MOSEK?) (x = [0 0 30 0 0], fval = -90)
clc
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
      1 0 1 -3 0];
b = [30;30];
lb = zeros(5,1);
ub = [40;1;inf;inf;1];
sos = '1';
sosind = [1:5]';
soswt = [1:5]';
sos.sostype = sos; sos.sosind = sosind; sos.soswt = soswt;

% [x,fval,exitflag,info] = opti_mintprog(f,A,b,[],[],lb,ub,[],sos,opts)

%% QP1 (x = [0.33 1.33 -0.66], fval = 2.833)
clc
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
rl = -Inf(3,1);
ru = [1;1;1];

[x,fval,exitflag,info] = mosekqp(H,f,A,rl,ru)

%% QP2 (x = [0.66 1.33], fval = -8.22)
clc
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1];
rl = -Inf(3,1);
ru = [2; 2; 3];
lb = [0;0];

[x,fval,exitflag,info] = mosekqp(H,f,A,rl,ru,lb)

%% QP3 (x = [0.3448 1.1034], fval = -6.4138)
clc
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2; 2 1; 1 1.5];
rl = [-Inf(3,1);2];
ru = [2; 2; 3; 2];
lb = [0;0];
ub = [10;10];

[x,fval,exitflag,info] = mosekqp(H,f,A,rl,ru,lb,ub)

%% QCQP1 (x = [0.0142 0.0048 0.1864], fval = -0.3776)
%note OPTI QC does not include 1/2 in front of x'Qx
clc
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
rl = -Inf(2,1);
ru = [20;30];
Q = eye(3);
l = [3;-2;5];
qrl = -Inf;
qru = 1;
lb = [0;0;0];
ub = [40;inf;inf];

[x,fval,exitflag,info] = mosekqcqp(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub)

%% QCQP2 (x = [0.5798 0.2899], fval = -1.7394)
%note OPTI QC does not include 1/2 in front of x'Qx
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
rl = -Inf(2,1);
ru = [2;5];
Q = {[1 0;0 1];[1 0;0 1]};
l = [0 2;2 -2];
qrl = -Inf(2,1);
qru = [1 1];
lb = [0;0];
ub = [40;inf];

opts = mosekset('display','iter');

[x,fval,exitflag,info] = mosekqcqp(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,[],opts)

%% MIQP1 (x = [1 2], fval = -11.5)
clc
%Objective & Constraints
H = [1 -1; -1 2];
f = -[2 6]';
A = [1 1; -1 2];
rl = -Inf(2,1);
ru = [3; 3.5];
lb = [0;0];
ub = [10;10];
xint = 'IC';

opts = mosekset('display','iter');

[x,fval,exitflag,info] = mosekmiqp(H,f,A,rl,ru,lb,ub,xint,[],opts) 

%% MIQP2 (x = [0.5 1 -0.5], fval = -2.75)
clc
%Objective & Constraints
H = eye(3);
f = -[2 3 1]';
A = [1 1 1;3 -2 -3; 1 -3 2]; 
rl = -Inf(3,1);
ru = [1;1;1];
xint = 'CIC';

[x,fval,exitflag,info] = mosekmiqp(H,f,A,rl,ru,[],[],xint)

%% MIQCQP1 (x = [0.0145 0.087 0], fval = -0.0942)
clc
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
rl = -Inf(2,1);
ru = [20;30];
Q = eye(3);
l = [0;0;0];
qrl = -Inf;
qru = 1;
lb = [0;0;0];
ub = [40;inf;inf];
xint = 'CCI';

mopt.MSK_IPAR_MIO_MAX_NUM_RELAXS = 10; %test field
opts = mosekset('display','iter','mskoption',mopt);

[x,fval,exitflag,info] = mosekmiqcqp(H,f,A,rl,ru,Q,l,qrl,qru,lb,ub,xint,[],opts)


%% SDP1
% http://docs.mosek.com/7.0/toolbox/A_guided_tour.html#section-node-_A%20guided%20tour_Semidefinite%20optimization
% https://www.controlengineering.co.nz/Wikis/OPTI/index.php/Probs/SDP
clc
clear prob
[r, res] = mosekopt('symbcon');

% min x1
% s.t. [x1 sqrt(2); sqrt(2) x1] >= 0

prob.c         = 0;

%C = [0 -sqrt(2); -sqrt(2) 0]
prob.bardim    = 2;
prob.barc.subj = 1; %which cone
prob.barc.subk = 2; %row
prob.barc.subl = 1; %col
prob.barc.val  = sqrt(2);

%A = [1 0; 0 1]
prob.a = sparse(1,1); 
prob.bara.subi = [1, 1]; %which A matrix
prob.bara.subj = [1, 1]; %which cone
prob.bara.subk = [1, 2]; %row
prob.bara.subl = [1, 2]; %col
prob.bara.val  = [1.0, 1.0];

prob.blc = 1;  %objective
prob.buc = 1; 

[r,res] = mosekopt('minimize info',prob); 

x = res.sol.itr.y
X = res.sol.itr.barx
S = res.sol.itr.bars

%% SDP2
clc
clear prob
[r, res] = mosekopt('symbcon');

% min x1+x2
% s.t. [x1 2; 2 x2] >= 0

prob.c         = [0;0];

%C = [0 2; 2 0]
prob.bardim    = 2;
prob.barc.subj = 1; %which cone
prob.barc.subk = 2; %row
prob.barc.subl = 1; %col
prob.barc.val  = 2;

%A1 = [1 0; 0 0]
%A2 = [0 0; 0 1]
prob.a = sparse(2,2); 
prob.bara.subi = [1, 2]; %which A matrix
prob.bara.subj = [1, 1]; %which cone
prob.bara.subk = [1, 2]; %row
prob.bara.subl = [1, 2]; %col
prob.bara.val  = [1.0, 1.0];

prob.blc = [1;1];  
prob.buc = [1;1]; 

prob.blx = [10;1];
prob.bux = [10;1];

[r,res] = mosekopt('minimize info',prob); 

x = res.sol.itr.y
X = res.sol.itr.barx
S = res.sol.itr.bars


%% SDP4
clc
clear prob
[r, res] = mosekopt('symbcon');

% min x1+x2+x3
% s.t. [x1 1 2; 1 x2 3; 2 3 100] >= 0

prob.c = [0;0;0];

%C = [0 1 2; 1 0 3; 2 3 100]
prob.bardim    = 3;
prob.barc.subj = [1 1 1 1]; %which cone
prob.barc.subk = [2 3 3 3]; %row
prob.barc.subl = [1 1 2 3]; %col
prob.barc.val  = [1 2 3 100];

%A1 = [1 0 0; 0 0 0; 0 0 0]
%A2 = [0 0 0; 0 1 0; 0 0 0]
%A3 = [0 0 0; 0 0 0; 0 0 0]
prob.a = sparse(3,3); 
prob.bara.subi = [1, 2]; %which A matrix
prob.bara.subj = [1, 1]; %which cone
prob.bara.subk = [1, 2]; %row
prob.bara.subl = [1, 2]; %col
prob.bara.val  = [1.0, 1.0];

prob.blc = [1;1;1];  
prob.buc = [1;1;1]; 

[r,res] = mosekopt('minimize info',prob); 

x = res.sol.itr.xx


%% Large Problems
clc
opts = optiset('solver','mosek','display','iter');
Opt = opti('control3.dat-s',opts)
[~,f,~,info] = solve(Opt)



