%% SOS Testing
%% LP Solve Example SOS1
clc
f = [-1 -1 -3 -2 -2]';
A = [-1 -1 1 1 0;
     1 0 1 -3 0];
b = [30;30];
ub = [40;1;inf;inf;1];

sostype = '12';
sosind = {[1:2]' [3:5]'};
soswt = {[1:2]' [3:5]'};

% [x,fval,ef,stat] = cplexmilp(f,A,b,[],[],sos,sosind,soswt,[],ub)

% SOS 1 OPTI Version
opts = optiset('solver','cbc');

prob = optiprob('f',f,'ineq',A,b,'ub',ub,'sos',sostype,sosind,soswt,'options',opts)

sos = [];
sos.type = sostype;
sos.index = sosind;
sos.weight = soswt;
prob2 = optiprob('f',f,'ineq',A,b,'ub',ub,'sos',sos,'options',opts)

Opt = opti(prob,opts)

[x,fval,ef,stat] = solve(Opt)

%% LP Solve Example SOS2
clc
f = [-1 -1 -3 -2 -2]';
A = sparse([-1 -1 1 1 0;
             1 0 1 -3 0]);
b = [30;30];

lb = zeros(5,1);
ub = [40;1;inf;inf;1];

sos = '1';
sosind = (1:5)';
soswt = (1:5)';

% [x,fval,ef,stat] = cplexmilp(f,A,b,[],[],sos,sosind,soswt,lb,ub,'IIIII');
% x
% stat

% SOS 2 OPTI Version
opts = optiset('solver','cbc');

prob = optiprob('f',f,'ineq',A,b,'bounds',lb,ub,'sos',sos,sosind,soswt,'options',opts)

Opt = opti(prob,opts)

[x,fval,ef,stat] = solve(Opt)


% [x,fval,ef,stat] = solve(Opt);
% x
% fval
% stat

%%
clc
% Objective (f'x)
f = [-1 -1 -3 -2 -2]';

% Row Constraints (rl <= A*x <= ru)
A = sparse([-1 -1  1  1  0;     %sparse A
             1  0  1 -3  0]);
rl = [-Inf;-Inf];    
ru = [30;30];

% Bounds (lb <= x <= ub)
lb = [0;0;0;0;0];
ub = [40;1;Inf;Inf;1];        

% SOS
sos_type = '1'; 
sos_index = [1 2 3 4 5]';
sos_weight = [1 2 3 4 5]';

% Setup Options
opts = optiset('solver','cbc'); %CBC is a SOS constraint solver

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'sos',sos_type,sos_index,...
           sos_weight,'options',opts)

% Solve Problem
[x,fval] = solve(Opt)

%%
% SOS Structure
sos = [];
sos.type = '12'; %each character represents a SOS set 
sos.index = {[1 2]' [3:5]'};
sos.weight = {[1 2]' [1:3]'};

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'sos',sos,'options',opts)

% Solve Problem
[x,fval] = solve(Opt)

