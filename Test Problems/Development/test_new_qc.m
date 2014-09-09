%% FIRST - OLD API

%% QCQP1 [-0.4004]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [20;30];
Q = eye(3);
l = [2;2;2];
r = 0.5;
lb = [0;0;0];
ub = [40;inf;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('display','iter','solver','ipopt'))
[x,fval,exitflag,info] = solve(Opt)

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('display','iter','solver','scip'))
[x,fval,exitflag,info] = solve(Opt)

%% QCQP3 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1]};
l = [0 2;2 -2];
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% QCQP3a [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
%Alternate QC form
Q = {[1 0;0 1];[1 0;0 1]};
l = {[0;2];[2 -2]};
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% QCQP4 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1];[1 0;0 1]};
l = [0 2 -2;2 -2 -1];
r = [1 1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% Positive Semidefinite QCQP
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;3.8];
Q = [0 0;0 10];
l = [0;-2];
r = 3;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('solver','auto'))
[x,fval,exitflag,info] = solve(Opt)

%% Indefinite QCQP
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [-1 0;0 1];
l = [0;-2];
r = -0.5;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qc',Q,l,r,'bounds',lb,ub,'options',optiset('solver','scip'))
[x,fval,exitflag,info] = solve(Opt)


%% NOW - NEW API

%% QCQP1 [-0.4004]
clc
%Objective & Constraints
H = [33 6    0;
     6  22   11.5;
     0  11.5 11];
f = [-1;-2;-3];
A = [-1 1 1; 1 -3 1];
b = [20;30];
Q = eye(3);
l = [2;2;2];
r = 0.5;
lb = [0;0;0];
ub = [40;inf;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf,r,'bounds',lb,ub,'options',optiset('display','iter','solver','scip'))
[x,fval,exitflag,info] = solve(Opt)

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf,r,'bounds',lb,ub,'options',optiset('display','iter','solver','scip'))
[x,fval,exitflag,info] = solve(Opt)

%% QCQP3 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1]};
l = [0 2;2 -2];
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf(2,1),r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% QCQP3a [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
%Alternate QC form
Q = {[1 0;0 1];[1 0;0 1]};
l = {[0;2];[2 -2]};
r = [1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf(2,1),r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% QCQP4 [-1.7394]
clc
H = zeros(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1];[1 0;0 1];[1 0;0 1]};
l = [0 2 -2;2 -2 -1];
r = [1 1 1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf(3,1),r,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

%% Positive Semidefinite QCQP
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;3.8];
Q = [0 0;0 10];
l = [0;-2];
r = 3;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,-Inf,r,'bounds',lb,ub,'options',optiset('solver','auto'))
[x,fval,exitflag,info] = solve(Opt)

%% Indefinite QCQP
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [-1 0;0 1];
l = [0;-2];
r = -0.5;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'Q',Q,'l',l,'qru',r,'bounds',lb,ub,'options',optiset('solver','scip'))
[x,fval,exitflag,info] = solve(Opt)

%% QUADRATIC EQUALITIES

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,r,r,'bounds',lb,ub,'options',optiset('display','iter','solver','auto'))
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)



%% QUADRATIC DOUBLE SIDED

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = [1 0;0 1];
l = [0;-2];
qrl = 3;
qru = 5;
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,qrl,qru,'bounds',lb,ub,'options',optiset('display','iter','solver','auto'))
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)

%% BOTH
%% QCQP3 [-2.8166]
clc
%Objective & Constraints
H = eye(2);
f = [-2;-2];
A = [-1 1; 1 3];
b = [2;5];
Q = {[1 0;0 1]; [1 0; 0 1]};
l = {[0;-2]; [-2;2]};
qrl = [3;1];
qru = [5;1];
lb = [0;0];
ub = [40;inf];
%Build & Solve
Opt = opti('H',H,'f',f,'ineq',A,b,'qcrow',Q,l,qrl,qru,'bounds',lb,ub,'options',optiset('display','iter','solver','ipopt'))
[x,fval,exitflag,info] = solve(Opt,[2;0])
plot(Opt)
