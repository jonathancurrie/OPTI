%% LIPSOL TESTING

%% LP1 [-31.4]
clc
%Objective & Constraints
f = -[6 5]';
A = ([1,4; 6,4; 2, -5]); 
b = [16;28;6];  
lb = [0;0];
ub = [10;10];
%Solve
[x,fval,ef,info] = opti_lipsol(f,A,b,[],[],lb,ub)

%% LP2
clc
%Objective & Constraints
f = -[1 2 3]';
A = [-1,1,1; 1,-3,1];
b = [20,30]';
Aeq = [1 1 1];
beq = 40;
lb = [0;0;0];
ub = [40;inf;inf];
%Solve
[x,fval,ef,info] = opti_lipsol(f,A,b,Aeq,beq,lb,ub)


%% LP3 (needs lower bounds...)
% clc
% %Objective & Constraints
% f = [8,1]';
% A = [-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]; 
% b = [-4,2,21,39,3,0,0]';
% %Solve
% [x,fval,ef,info] = opti_lipsol(f,A,b)

%% Raw LIPSOL Test
clc
%Problem
A = sparse([ 1     3     4     1     0
             5    -1     1     0    -1
             2     1    -1     0     0]);
b = [8;0;4];
c = [2;1;-2.5;0;0];
lb = [2;0;0;0;0];
ub = [1e32;1e32;9;1e32;1e32];
[xsol,fp,fd,info,msg,times] = lipsol(A,b,c,lb,ub)

%% MPS Files
clc
prob = coinRead('maros-r7.mps');
opts = optiset('solver','lipsol','display','iter','maxiter',1e4);
Opt = opti(prob,opts);
%Solve
[~,f,e,i] = solve(Opt)




