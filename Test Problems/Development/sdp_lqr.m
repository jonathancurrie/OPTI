%% Solving LQR via SDP (LMIs)
% http://www1.se.cuhk.edu.hk/~zhang/Courses/Seg5660/Lecture%206.pdf

%% Continuous LQR
clc
Gs = tf(2,[0.5 0.2 0.1]);
[A,B] = ssdata(ss(Gs));
Q = diag([10;10]); R = 1;

[K,S] = lqr(A,B,Q,R)

%LQR LMI - Vars are P11 P12 P22
P = sym('[P11 P12; P12 P22]')
Lp = vpa([R B'*P; P*B Q + A'*P + P*A],5)

% Objective (max)
f = -[1 0 1]; %just diagonal elements

% Semidefinite Constraint (cheap way to pick off elements)
C = -double(subs(Lp,P,zeros(2)));
A1 = double(diff(Lp,P(1)));
A2 = double(diff(Lp,P(2)));
A3 = double(diff(Lp,P(end)));
sdcone = sparse([C(:) A1(:) A2(:) A3(:)]);

% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone)
% Solve the SDP
x = solve(Opt)

%Reform P and solve K
P = [x(1) x(2); x(2) x(3)];
Ksdp = R\(B'*P) %see help lqr

%% Discrete LQR
clc
Gd = c2d(Gs,0.1);
[A,B] = ssdata(ss(Gd));
Q = diag([10;10]); R = 1;

[K,S] = dlqr(A,B,Q,R)

%LQR LMI - Vars are P11 P12 P22
P = sym('P',[2 2]); 
Lp = vpa([B'*P*B+R B'*P*A; A'*P*B Q + A'*P*A - P],5)

% Objective (max)
f = -[1 0 1]; %just diagonal elements

% Semidefinite Constraint (cheap way to pick off elements)
C = -double(subs(Lp,P,zeros(2)));
A1 = double(diff(Lp,P(1)));
A2 = double(diff(Lp,P(2))) + double(diff(Lp,P(3)));
A3 = double(diff(Lp,P(end)));

sdcone = sparse([C(:) A1(:) A2(:) A3(:)]);

% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone)
% Solve the SDP
x = solve(Opt)

%Reform P and solve K
P = [x(1) x(2); x(2) x(3)];
Ksdp = (B'*P*B + R)\(B'*P*A) %see help dlqr

%% Continuous Lyapunov
clc
[A,B] = ssdata(ss(Gs));
X = lyap(A,Q)

%LYAP LMI - Vars are P11 P12 P22
P = sym('P',[2 2]); 
Lp = vpa(A*P + P*A' + Q)

% Objective (max)
f = -[1 1 1]; %just diagonal elements

% Semidefinite Constraint (cheap way to pick off elements)
C = -double(subs(Lp,P,zeros(2)));
A1 = double(diff(Lp,P(1)));
A2 = double(diff(Lp,P(2))) + double(diff(Lp,P(3)));
A3 = double(diff(Lp,P(end)));

sdcone = sparse([C(:) A1(:) A2(:) A3(:)]);

% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone,'solver','dsdp')
% Solve the SDP
[x,f,e,i] = solve(Opt)
P = [x(1) x(2); x(2) x(3)];

A*P + P*A' + Q
A*X + X*A' + Q

%% Discrete Lyapunov
clc
[A,B] = ssdata(ss(Gd));
X = dlyap(A,Q)

%LYAP LMI - Vars are P11 P12 P22
P = sym('P',[2 2]); 
Lp = vpa(A*P*A' - P + Q)

% Objective (max)
f = -[1 1 1]; %just diagonal elements

% Semidefinite Constraint (cheap way to pick off elements)
C = -double(subs(Lp,P,zeros(2)));
A1 = double(diff(Lp,P(1)));
A2 = double(diff(Lp,P(2))) + double(diff(Lp,P(3)));
A3 = double(diff(Lp,P(end)));

sdcone = sparse([C(:) A1(:) A2(:) A3(:)]);

% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone,'solver','dsdp')
% Solve the SDP
[x,f,e,i] = solve(Opt)
P = [x(1) x(2); x(2) x(3)];

A*P*A' - P + Q
A*X*A' - X + Q

