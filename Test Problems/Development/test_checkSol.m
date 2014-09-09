%% Checking Solutions
% note remove the private setting from the opti class
clear all

%% Bounds
clc
f = [-1;-2;-3;-4];
lb = [0;0;0;0];
ub = [1;1;1;1];

O = opti('f',f,'bounds',lb,ub,'options',optiset('solver','clp'));
O.sol = [-1;-1;2;2];
O.ef = 1;

[status,msg] = checkSol(O)

%% LP General Constraints
clc
f = [-1;-2];
A = [1 0;0 1];
b = [0;0];
Aeq = [1 0; 0 1];
beq = [0.5;0.5];

O = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'options',optiset('solver','cplex'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% LP Row Constraints
clc
f = [-1;-2];
A = sparse([1 0;0 1;1 0;0 1;1 0;0 1]);
rl = [2;2;-Inf;-Inf;0.5;0.5];
ru = [Inf;Inf;0;0;0.5;0.5];

O = opti('f',f,'lin',A,rl,ru,'options',optiset('solver','clp'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% QC Single Constraint
clc
f = [-1;-2];
Q = eye(2);
l = [-1;-2];

O = opti('H',zeros(2),'f',f,'qc',Q,l,-1.5,'options',optiset('solver','scip'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% QC Multiple Constraints
clc
f = [-1;-2];
Q = {eye(2);eye(2)};
l = [-1 -1;-2 -3];

O = opti('H',zeros(2),'f',f,'qc',Q,l,[-1.5;-3],'options',optiset('solver','scip'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% NLP General Constraints
clc
fun = @(x) x(1)+x(2);
nlcon = @(x) [x(1);x(2);x(1);x(2);x(1);x(2)];
nlrhs = [2;2;0;0;0.5;0.5];
nle = [1;1;-1;-1;0;0];

O = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ndec',2,'options',optiset('solver','nlopt'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% NLP Row Constraints
clc
fun = @(x) x(1)+x(2);
nlcon = @(x) [x(1);x(2);x(1);x(2);x(1);x(2)];
cl = [2;2;-Inf;-Inf;0.5;0.5];
cu = [Inf;Inf;0;0;0.5;0.5];

O = opti('fun',fun,'nl',nlcon,cl,cu,'ndec',2,'options',optiset('solver','ipopt'));
O.sol = [1;1];
O.ef = 1;

[status,msg] = checkSol(O)

%% Integer Constraints
clc
f = [-1;-2;-3;-4];
lb = [0;0;0;0];
ub = [1;1;1;1];

O = opti('f',f,'bounds',lb,ub,'int','ICCI','options',optiset('solver','cbc'));
O.sol = [-1.1;-1.5;2;2.1];
O.ef = 1;

[status,msg] = checkSol(O)

%% Binary Constraints
clc
f = [-1;-2;-3;-4];
lb = [0;0;0;0];
ub = [1;1;1;1];

O = opti('f',f,'bounds',lb,ub,'int','ICBI','options',optiset('solver','cbc'));
O.sol = [-1.1;-1.5;1.001;2.1];
O.ef = 1;

[status,msg] = checkSol(O)