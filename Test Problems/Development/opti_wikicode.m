%% OPTI Wiki Code
% A copy of all code on the Wiki that is run before each release to ensure
% compatibility. Also regenerates all plots if required.
clc
clear
clf
%Setup Figure Generation
% Remember to crop all spy plots (3x) + higher D plots (2x)!
figDir = ''; %'C:\OPTIFigs'; %empty to not generate plots
figure(1);
stdsize = [680   558   560   420];
set(gcf,'position',stdsize);
set(gcf,'PaperPositionMode','auto');

%% Basic OPTI Usage
% https://www.inverseproblem.co.nz/OPTI/index.php/GetStart/Basics

% Ex1
% Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Create OPTI Object
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub)
[x,fval,exitflag,info] = solve(Opt)

% Alt Strategies
Opt = opti('c',f,'ineq',A,b,'bounds',lb,ub);    %c = f

% OR
Opt = opti('grad',f,'ineq',A,b,'bounds',lb,ub); %grad = f

% OR
Opt = opti('f',f,'A',A,'b',b,'bounds',lb,ub);   %individual A,b

% OR
Opt = opti('f',f,'ineq',A,b,'lb',lb,'ub',ub);   %individual lb,ub

% Choosing a Solver
optiSolver
optiSolver('matrix')
optiSolver('LP')
optiSolver('CLP')
optiset

%Specify to use the COIN-OR LP solver CLP
opts = optiset('solver','clp');

%Rebuild the problem passing the options via the 'options' argument
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'options',opts)

optiSolver('config')

% Ex2
% Objective
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2; 

% Constraints
lb = [-5;-5];
ub = [5;5];

% Initial Starting Guess (Required for Nonlinear Problems)
x0 = [0;0];

% Setup Options
opts = optiset('solver','lbfgsb','display','iter');

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)

% Solve!
[x,fval,ef,info] = solve(Opt,x0)
plot(Opt,7,1)  %plot the object, +-7 from the solution and take log(obj)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex2_plot']);
end

%% MATLAB Optimization Toolbox Overloads
% https://www.inverseproblem.co.nz/OPTI/index.php/GetStart/Overloads

% Ex1
% Problem
f = -[6 5]';                %Objective Function (min f'x)
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Solve using MATLAB's linprog:
if (exist('linprog.m','file'))
    x = linprog(f,A,b,[],[],lb,ub)  
end

% Solve using corresponding OPTI overload:
x = opti_linprog(f,A,b,[],[],lb,ub)

% Ex2
% Objective
f = [2, 3, 7, 7]';         
% Constraints       
A = [-1, -1, 2, 5;1, -2, -1, -4];
b = [-2; -3];  
lb = zeros(4,1); 
ub = [30 100 20 1]';  

% Integer Constraints
xtype = [2 4];  %variables 1, 3 are continuous, 2, 4 are integer

% Solve using OPTI overload:
x = opti_intlinprog(f,xtype,A,b,[],[],lb,ub)

% Ex3 NLP
% Objective
fun = @(x) sin(x(1) + x(2)) + (x(1) - x(2))^2 - 1.5*x(1) + 2.5*x(2) + 1;          
% Constraints       
lb = [-1.5;-3];
ub = [4;3];  

% Starting Guess
x0 = [0;0];

% Setup Options
opts = optiset('solver','nomad');

% Solve using OPTI overload:
[x,fval] = opti_fmincon(fun,x0,[],[],[],[],lb,ub,[],opts)

%% Linear Program
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/LP

% Ex 1
% Objective
f = -[6 5]';                %Objective Function Vector (min f'x)

% Constraints
A = [1,4; 6,4; 2,-5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Create OPTI Object
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub)

% Solve the LP problem
[x,fval,exitflag,info] = solve(Opt)

plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exlp']);
end

% Ex2 
% Load the LP from .mps file
prob = coinRead('maros-r7.mps')

% Examine A matrix
spy(prob.A)

if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exlpspy']);
end

prob

% Create an OPTI object from the problem structure
Opt = opti(prob)

% Solve the resulting model
[x,fval] = solve(Opt);

% Ex3
% Objective (f'x)
f = -[1 2 3]';

% Row Constraints (rl <= A*x <= ru)
A = sparse([-1  1  1;     %sparse even though all nz
             1 -3  1;
             1  1  1]);
rl = [-Inf;-Inf;40];      %top two rows are only Ax <= b
ru = [20;30;40];

% Bounds (lb <= x <= ub)
lb = [0;0;0];
ub = [40;Inf;Inf];        %x2 and x3 are unbounded above

% Setup Options
opts = optiset('solver','clp'); %CLP is a row constraint solver

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'options',opts)

% Solve Problem
[x,fval] = solve(Opt)

% Ex4
% Objective (f'x)
f = -[1 2 3]';

% Linear Constraints
A = sparse([-1  1  1;     %sparse even though all nz
            -1  3 -1;     %note this row inverted just for example purposes
             1  1  1]);
b = [20;-30;40];
e = [-1;1;0]; %-1 <=, 0 ==, 1 >=

% Bounds (lb <= x <= ub)
lb = [0;0;0];
ub = [40;Inf;Inf];        %x2 and x3 are unbounded above

% Setup Options
opts = optiset('solver','clp'); %CLP is a row constraint solver

% Build OPTI Object
Opt = opti('f',f,'mix',A,b,e,'bounds',lb,ub,'options',opts)

% Solve Problem
[x,fval] = solve(Opt)

% Ex5
% Objective (f'x + objbias)
f = -[6 5]';
objbias = 5;

% Linear Constraints (A*x <= b)
A = [1,4; 6,4; 2, -5]; 
b = [16;28;6];    

% Bounds (lb <= x <= ub)
lb = [0;0]; ub = [10;10];

% Build OPTI Object
Opt = opti('f',f,'objbias',objbias,'ineq',A,b,'bounds',lb,ub)

% Solve Problem
[x,fval] = solve(Opt)

%% MILP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/MILP

% Ex1
% Objective
f = -[6 5]';                %Objective Function Vector (min f'x)

% Constraints
A = [1,4; 6,4; 2, -5];      %Linear Inequality Constraints (Ax <= b)
b = [16;28;6];    
lb = [0;0];                 %Bounds on x (lb <= x <= ub)
ub = [10;10];

% Integer Constraints
xtype = 'II';               %x1 & x2 are Integer
% Create OPTI Object
Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'xtype',xtype)

% Solve the MILP problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exmilp']);
end

%Ex 2
xtype = 'CBI';
xtype = [2 3];

% Ex3
% Load the MILP from .mps file
prob = coinRead('testMILP2.mps')

% Examine A matrix
spy(prob.A)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exmilpspy']);
end

% Create an OPTI object from the problem structure
Opt = opti(prob)

% Solve the resulting model
[x,fval] = solve(Opt)

% Ex4
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


% SOS Structure
sos.type = '12'; %each character represents a SOS set 
sos.index = {[1 2]' [3:5]'};
sos.weight = {[1 2]' [1:3]'};

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub,'sos',sos,'options',opts)

% Ex5
% Objective
nC = 10; %Number of Continuous Variables
nI = 10; %Number of Integer Variables
nB = 10; %Number of Binary Variables

% Build xtype vector
xtype = [repmat('C',1,nC),repmat('I',1,nI),repmat('B',1,nB)]

%% SDP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/SDP

% Ex1
% Objective
f = 1;                

% Semidefinite Constraint
F0 = -[0 sqrt(2); sqrt(2) 0];
F1 = eye(2);
sdcone = sparse([F0(:) F1(:)]);
% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone)

% Solve the QP problem
[x,fval,exitflag,info] = solve(Opt)

% Ex2
% Objective
f = [1;1];
% Bounds on x (lb <= x <= ub)
lb = [0;0];            
ub = [10;10];
% Semidefinite Constraint (using alternative notation)
C = -[0 2; 2 0];
A1 = [1 0; 0 0];
A2 = [0 0; 0 1];
sdcone = sparse([C(:) A1(:) A2(:)]);

% Options
opts = optiset('solver','dsdp','display','iter');

% Create OPTI Object
Opt = opti('f',f,'bounds',lb,ub,'sdcone',sdcone,'options',opts)

% Solve the bounded SDP
[x,f] = solve(Opt)
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2sdp']);
end

% Ex 3
clear sdcone
% Objective
f = [1;0;0;0];
% Semidefinite Constraint 1
F0 = zeros(2);
F1 = eye(2);
F2 = -[1 0; 0 0];
F3 = -[0 1; 1 0];
F4 = -[0 0; 0 1];
sdcone{1} = [F0 F1 F2 F3 F4];
% Semidefinite Constraint 2
F0 = [1 0.2; 0.2 1];
F1 = zeros(2);
F2 = [1 0; 0 0];
F3 = [0 1; 1 0];
F4 = [0 0; 0 1];
sdcone{2} = [F0 F1 F2 F3 F4];

% Options
opts = optiset('solver','csdp','display','iter');

% Create OPTI Object
Opt = opti('f',f,'sdcone',sdcone,'options',opts)

% Solve the SDP
[x,f] = solve(Opt)

%Ex 4
% Load SeDuMi problem variables from MAT file
load sdp_truss1.mat

% Options
opts = optiset('display','iter');

% Create OPTI Object from SeDuMi Args
Opt = opti('sedumi',At,b,c,K,'options',opts)

% Solve the SDP
if(~isempty(which('sedumi.m')))
    [x,f] = solve(Opt)
end

%% QP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/QP

% Ex 1
% Objective
H = [1 -1; -1  2];          %Objective Function (min 0.5x'Hx + f'x)
f = -[2 6]';                

% Constraints
A = [1,1; -1,2; 2, 1];      %Linear Inequality Constraints (Ax <= b)
b = [2;2;3];    
lb = [0;0];                 %Bounds on x (lb <= x)
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb)

% Solve the QP problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exqp']);
end

%Ex2
% Load the QP from .qps file
prob = coinRead('MPCqp1.qps')

% Examine A matrix
spy(prob.A')
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2qpspy']);
end

% Create an OPTI object from the problem structure
Opt = opti(prob)

% Solve the resulting model
[x,fval] = solve(Opt)


e = eig(H)

% Ex 3 Plots
data.d = [1 1]; data.name = 'Positive Definite QP';
opti_wikiplot('qp',data)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_pdqp']);
end

data.d = [-1 -1]; data.name = 'Negative Definite QP';
opti_wikiplot('qp',data)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ndqp']);
end

data.d = [0 1]; data.name = 'Positive Semi-Definite QP';
opti_wikiplot('qp',data)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_psdqp']);
end

data.d = [-1 1]; data.name = 'Indefinite QP';
opti_wikiplot('qp',data)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_idqp']);
end

% Ex 4
%Problem
H = [0 -2; -2 0];           %Objective Function (min 1/2x'Hx + f'x)
f = [0 0]';
lb = [-0.5;-0.5];           %Bounds on x (lb <= x <= ub)   
ub = [1;1];

% Options
opts = optiset('solver','clp');

% Create OPTI Object
Opt = opti('qp',H,f,'bounds',lb,ub,'options',opts)

% Solve the Indefinite QP
[x,f,~,info] = solve(Opt)

plot(Opt,1.5)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2qp']);
end

%% MIQP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/MIQP

% Ex 1
% Objective
H = [1 -1; -1  2];          %Objective Function (min 0.5x'Hx + f'x)
f = -[2 6]';                

% Constraints
A = [1,1; -1,2; 2, 1];      %Linear Inequality Constraints (Ax <= b)
b = [2;2;3];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Integer Constraints
xtype = 'IC';

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'xtype',xtype)

% Solve the MIQP problem
[x,fval,exitflag,info] = solve(Opt)

plot(Opt,2)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exmiqp']);
end

%% QCQP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/QCQP

% Ex1
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraint
Q = [1 0; 0 1];             %Quadratic Inequality (x'Qx + l'x <= r)
l = [0;-2];
r = 1;
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exqcqp']);
end

%Ex 2
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Inequalities (x'Qx + l'x <= r)
     [1 0; 0 1]};
l = [[0;-2] [-1;2]];
r = [1;1.2];

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2qcqp']);
end

%Ex3
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Inequalities (x'Qx + l'x <= r)
     [1 0; 0 1]};
l = {[0;-2]; [-1;2]};
r = {1; 1.2};

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)

% Ex4
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;3.8];    
lb = [0;0];                 %Bounds on x (lb <= x)
ub = [40;inf];

% Quadratic Constraints
Q = [0 0; 0 10];            %Quadratic Inequalities (x'Qx + l'x <= r)
l = [0;-2];
r = 3

% Set OPTI Options
opts = optiset('solver','ipopt','display','iter');

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'bounds',lb,ub,'qc',Q,l,r,'options',opts)

% Solve the QCQP problem [as an NLP]
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt,3)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex4qcqp']);
end

% Ex5
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Linear Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraints
Q = {[1 0; 0 1]             %Quadratic Constraints (qrl <= x'Qx + l'x <= qru)
     [1 0; 0 1]};
l = {[0;-2]; [-2;2]};
qrl = {3; 1};               %QC1 is double sided, QC2 is an equality
qru = {5; 1};

% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qcrow',Q,l,qrl,qru)

% Solve the QCQP problem
[x,fval,exitflag,info] = solve(Opt)

% Plot Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex5qcqp']);
end

%% MIQCQPs
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/MIQCQP

%Ex1
% Objective
H = eye(2);                 %Objective Function (min 0.5x'Hx + f'x)
f = -[2 2]';                

% Constraints
A = [-1,1; 1,3];            %Linear Inequality Constraints (Ax <= b)
b = [2;5];    
lb = [0;0];                 %Bounds on x (lb <= x)

% Quadratic Constraint
Q = [1 0; 0 1];             %Quadratic Inequality (x'Qx + l'x <= r)
l = [0;-2];
r = 1;

% Integer Constraints
xtype = 'IC';
% Create OPTI Object
Opt = opti('qp',H,f,'ineq',A,b,'lb',lb,'qc',Q,l,r,'xtype',xtype)

% Solve the MIQCQP problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,3)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exmiqcqp']);
end

%% SNLE
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/SNLE

% Ex 1
% System of Nonlinear Equations
nleq = @(x) [ 2*x(1) - x(2) - exp(-x(1));
            -x(1) + 2*x(2) - exp(-x(2))];

% Starting Guess
x0 = [-5;5];
% Create OPTI Object
Opt = opti('nleq',nleq,'x0',x0)

% Solve the SNLE problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,[],1)

optiSolver('SNLE')

% Ex2
% System of Nonlinear Equations
nleq = @(x) [10*(x(2) - x(1)^2)
            sqrt(90)*(x(4) - x(3)^2)
            sqrt(10)*(x(2) + x(4) - 2)
            (1/sqrt(10))*(x(2) - x(4))];

% Nonlinear Equations Jacobian
nlJac = @(x) sparse([-20*x(1),10,0,0
                     0,0,-6*10^(1/2)*x(3),3*10^(1/2)
                     0,10^(1/2),0,10^(1/2)
                     0,10^(1/2)/10,0,-10^(1/2)/10]);

% Jacobian Sparsity Pattern
nlJacstr = @() sparse([1 1 0 0
                       0 0 1 1
                       0 1 0 1
                       0 1 0 1]);

% Starting Guess
x0 = [-30;-10;-30;-10];

% Sparse SNLE OPTI Problem
Opt = opti('nleq',nleq,'nlJac',nlJac,'nlJacstr',nlJacstr,'x0',x0)

% Solve
[x,fval,exitflag,info] = solve(Opt)

%% SCNLE
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/SCNLE

% Ex1
% System of Nonlinear Equations
nleq = @(x) [ 2*x(1) - x(2) - exp(-x(1));
            -x(1) + 2*x(2) - exp(-x(2))];

% Bounds
lb = [0.6;0];
ub = [1;1];

% Starting Guess
x0 = [-5;5];
% Create OPTI Object
Opt = opti('nleq',nleq,'bounds',lb,ub,'x0',x0)

% Solve the SCNLE problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,0.75,1)

%% NLS
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/NLS

% Ex1
% Objective (Fitting) Function
fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];

% Fitting Data
ydata = [0;0];

% Starting Guess
x0 = [-1.2;1];
% Create OPTI Object
Opt = opti('fun',fun,'ydata',ydata,'x0',x0)

% Solve the NLS problem
[x,fval,exitflag,info] = solve(Opt)
plot(Opt,[],1) %tell OPTI to take the log of the objective for detail

% Ex2
% Objective (Fitting) Function
fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];

% Fitting Data
ydata = [0;0];

% Bounds
lb = [-2;-2];
ub = [0.5;0.5];

% Starting Guess
x0 = [-1.2;-2];

% Create OPTI Object
Opt = opti('fun',fun,'ydata',ydata,'bounds',lb,ub,'x0',x0)

% Solve using MKLTRNLS (auto-selected)
[x,fval,exitflag,info] = solve(Opt)

% Ex3
% Objective (Fitting) Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);

% Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]';

% Setup Options
opts = optiset('solver','lmder','display','iter');

% Starting Guess
x0 = [100; -1];

% Create OPTI Object
Opt = opti('fun',fun,'data',xdata,ydata,'x0',x0,'options',opts)

% Solve
[x,fval,exitflag,info] = solve(Opt)
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex3nls']);
end

% Ex 4
% Objective (Fitting) Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);

% Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742,...
       3.0237, 2.7002, 2.8781, 2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271,...
       1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,1.3825, 1.5087, 1.3624,...
       1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];;

% Weighting Vector
wts = ones(size(ydata)); wts(end) = 1e3;

% Starting Guess
x0 = [1.0; 0.0; 0.0];

% Create OPTI Objects for both Weighted and Un-weighted Cases
OptW = opti('fun',fun,'ydata',ydata,'weights',wts,'x0',x0);
Opt = opti('fun',fun,'ydata',ydata,'x0',x0);

% Solve Each Case
xw = solve(OptW)
x = solve(Opt)

% Plot Comparison
plot(i,ydata,'ko',i(end),ydata(end),'ksq',i,fun(xw),'r*-',i,fun(x),'m*-')
xlim([length(i)*0.45 length(i)*1.01]); xlabel('i'); ylabel('y');
legend('Original Data','Target Point','Weighted Fit','Standard Fit');
title('NLS Curve Fit - Comparison of Weighted and Un-Weighted');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex4nlswt']);
end

%Ex 5
% Objective (Fitting) Function
fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];

% Objective Gradient (Matrix in NLS case)
grad = @(x) [-200*x(1),100; -1,0];

% Fitting Data
ydata = [0;0];

% Starting Guess
x0 = [-1.2;1];
% Create OPTI Object
Opt = opti('fun',fun,'grad',grad,'ydata',ydata,'x0',x0)

% Solve the NLS problem
[x,fval,exitflag,info] = solve(Opt)

% Ex6
% Fitting Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);

% Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]';

% Setup Options
opts = optiset('solver','nomad','display','iter');

% Starting Guess
x0 = [100; -1];

% Create OPTI Object
Opt = opti('fun',fun,'data',xdata,ydata,'x0',x0,'options',opts)

% Solve
[x,fval,exitflag,info] = solve(Opt)

%Ex7
% Fitting Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); 

% Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate

% Starting Guess
x0 = [5e-3; 2e-2];

% Create OPTI Object
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0)

% Solve
[x,fval,exitflag,info] = solve(Opt)

% Calculate fit statistics
fitStats(Opt)

% Plot solution with statistics
plot(Opt)
h = legend;
set(h,'location','northwest');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_fitstats']);
end

%% DNLS
% https://www.inverseproblem.co.nz/OPTI/index.php/Dynamic/DynamicSystemParameterEstimation

%Ex1
% ODE System
ode = @(t,z,p) [-p(1)*z(1) + 4; 
                2*z(1) - p(1)*z(2) + 5; 
                -4*z(1) - 2*z(2)*z(3) - p(2)];

% Initial Conditions
z0 = [-1.5;1.25;1];

% True Parameter Values
p = [2.345;1.1];

% Generate Fitting Data
tm  = 0:0.1:2;                 %measurement times
odeInt = @(t,z) ode(t,z,p);    %ODE function for ODE45
[~,zm] = ode45(odeInt,tm,z0);  %Solve ODEs
% Starting Guess
theta0 = [1;0.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex1dnls']);
end

%Ex2
% Indicate Initial Conditions to Solve for Using NaNs
z0 = [-1.5;NaN;NaN];
% Starting Guess [p;z0]
theta0 = [1;0.5;0.5;0.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2dnls']);
end

%Ex3
% Add Bounds
lb = [1.5;0;0.1;0.1];
ub = [3.5;2;1.5;1.5];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'bounds',lb,ub,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

%Ex4
% Set Dynamic Options (specifying states of interest)
dopts = optidynset('stateIndex',[2 3]);

% Add to General OPTI Options
opts = optiset('display','iter','dynamicOpts',dopts);

% Index Measured Data
zm23 = zm(:,[2 3]);
% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm23,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex3dnls']);
end

%Ex5
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm2  = 0:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
[~,zm2] = ode45(odeInt,tm2,z0); zm2 = zm2(:,2);
% Note 0 is required below just to match initial condition in this example
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3(2:end,3); %drop first point

% Group measurements and time stamps in cell arrays
tm_multi = {tm2;tm3};
zm_multi = {zm2;zm3};
% Return z0 to include NaNs of estimated initial conditions
z0 = [-1.5;NaN;NaN];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex4dnls']);
end

%Ex6
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm1  = 0.5:0.1:2;              %state 1
tm2  = 1:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
tm_multi = {tm1;tm2;tm3};
zm_multi = {zm1;zm2;zm3};
% We will need to estimate all initial conditions in this problem
z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
theta0 = [1;0.5;0.5;0.5;0.5];

% Tell OPTI we want to solve from t = 0
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex6tdnls']);
end

%Ex7
% Initial Conditions (re-initialize)
z0 = [-1.5;1.25;1];

% New ODE with modified parameters (representing a different run)
odeIntB = @(t,z) ode(t,z,p*0.85); 

% Measurement Times for each State
tm1  = 0.5:0.1:2;              %state 1
tm1b = 1:0.25:2;               %state 1 2nd run
tm2  = 1:0.1:2;                %state 2
tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1]  = ode45(odeInt, [0 tm1], z0);  zm1  = zm1((2:end),1);
[~,zm1b] = ode45(odeIntB,[0 tm1b],z0);  zm1b = zm1b((2:end),1); %2nd run
[~,zm2]  = ode45(odeInt, [0 tm2], z0);  zm2  = zm2((2:end),2);
[~,zm3]  = ode45(odeInt, [0 tm3], z0);  zm3  = zm3((2:end),3);
[~,zm3b] = ode45(odeInt, [0 tm3], z0*0.75); zm3b = zm3b((2:end),3); %2nd run

% Group measurements and time stamps in cell arrays
tm_multi = {[tm1 tm1b];tm2;[tm3 tm3]}; %concatenate repeated points
zm_multi = {[zm1;zm1b];zm2;[zm3;zm3b]};
% We will need to estimate all initial conditions in this problem
z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
theta0 = [1;0.5;0.5;0.5;0.5];

% Tell OPTI we want to solve from t = 0
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
           'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex7dnls']);
end

%Ex8
if (exist('syms.m'))
    % Declare Symbolic Variables
    syms p1 p2 z1 z2 z3

    % Symbolic ODE RHS
    ODE = [-p1*z1 + 4; 
           2*z1 - p1*z2 + 5; 
           -4*z1 - 2*z2*z3 - p2];

    % Solve Jacobians
    dfdz_sym = jacobian(ODE,[z1 z2 z3])
    dfdp_sym = jacobian(ODE,[p1 p2])
    % ODE System
     ode = @(t,z,p) [-p(1)*z(1) + 4; 
                     2*z(1) - p(1)*z(2) + 5; 
                     -4*z(1) - 2*z(2)*z(3) - p(2)];

    % Simple Symbolic + String Parsing Derivative Generation
    [dfdz,dfdp] = symDynJac(ode)
    % Analytical Derivative Expressions
    dfdz = @(t,z,p) [-p(1) 0,       0;
                     2,    -p(1),   0;
                     -4,   -2*z(3), -2*z(2)];
    dfdp = @(t,z,p) [-z(1), 0;
                     -z(2), 0;
                     0,    -1];

    % Set Dynamic Options
    dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0);

    % General OPTI Options
    opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);

    % Create OPTI Object
     Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
                'options',opts)

    % Solve
    [theta,fval,exitflag,info] = solve(Opt)

    % Plot the Solution
    plot(Opt)
end

% Ex9
% Flame ODE System
ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;

% Initial Condition (Controls Stiffness)
z0 = 0.001;

% True Parameter
p = 1.5;

% Generate Fitting Data
tm  = [0 65:70 80:0.1:90 100]*10;  %measurement times
odeInt = @(t,z) ode(t,z,p);        %ODE function for ODE15S
[~,zm] = ode15s(odeInt,tm,z0);     %Solve ODE
% Starting Guess
theta0 = 1;

% Set Options (use ode15s)
dopts = optidynset('sensitivity','none');
opts1 = optiset('display','iter','iterfun',@optiplotlogfval,'dynamicOpts',dopts);
opts2 = optiset('display','iter','iterfun',@optiplotlogfval);

% Create OPTI Objects
OptNoDer = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts1)
OptWDer = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts2)

% Solve
[theta,fval,exitflag,info] = solve(OptNoDer)
[theta,fval,exitflag,info] = solve(OptWDer)

% Plot the Solution
subplot(121); plot(OptNoDer); title('No Sensitivity');
subplot(122); plot(OptWDer); title('With Sensitivity');
clf %NOT WORKING IN 2014B - well, working too well...

clc
% Partial Derivatives
dfdz = @(t,z,p) 2*p(1)*z(1) - 3*z(1)^2;
dfdp = @(t,z,p) z(1)^2;

% Set Dynamic Options
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'integrator','ode15s',...
                   'sensitivity','User');

% General OPTI Options
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)


% Set Dynamic Options
dopts = optidynset('integrator','ode15s','sensitivity','None');

% General OPTI Options
opts = optiset('solver','nomad','display','iter','dynamicOpts',dopts);

% Create OPTI Object (always set finite bounds for derivative free)
Opt = opti('ode',ode,'data',tm,zm,'bounds',0.5,2.5,'z0',z0,...
           'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

%% NLP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/NLP

%Ex1
% Objective
fun = @(x) log(1 + x(1)^2) - x(2);    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2; %Nonlinear Equality Constraint
nlrhs = 4;
nle = 0;                              %Constraint type: -1 <=, 0 ==, 1 >=

% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ndec',2)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)

plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_exnlp']);
end

%Ex2
% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Nonlinear Constraints (cl <= nlcon(x) <= cu)
nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
               x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
cl = [25;40];
cu = [Inf;40];

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

% Initial Guess
x0 = [1 5 5 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

%Ex3
% Objective
fun = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;

% Linear Equality Constraints
Aeq = [1 3 0 0 0;
       0 0 1 1 -2;
       0 1 0 0 -1];
beq = [4;0;0];

% Starting Guess
x0 = [ 2.5 0.5 2 -1 0.5 ];

% Options
opts = optiset('solver','ipopt');

% Create OPTI Object
Opt = opti('fun',fun,'eq',Aeq,beq,'options',opts)

% Solve using IPOPT
[x,fval,exitflag,info] = solve(Opt,x0)

%Ex4
% Objective
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

% Starting Guess
x0 = [0;0];

% Create OPTI Object
Opt = opti('fun',fun,'x0',x0)

% Solve
[x,fval,exitflag,info] = solve(Opt)

%Ex5
% Objective
fun = @(x) log(1 + x(1)^2) - x(2);

% Objective Gradient (row vector)
grad = @(x) [2*x(1)/(x(1)^2+1) -1];

% Nonlinear Constraints
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2; 
nlrhs = 4;
nle = 0;       

% Constraint Jacobian (matrix, one row per constraint)
jac = @(x) [4*x(1)*(x(1)^2+1), 2*x(2)];
% Create OPTI Object
Opt = opti('fun',fun,'grad',grad,'nlmix',nlcon,nlrhs,nle,'jac',jac,'ndec',2)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)

% Ex6
% Objective
fun = @RiemannND;

% Bounds
lb = 0.5;
ub = 2;

% Starting Guess
x0 = 0.5;

% Options
opts = optiset('solver','ipopt');

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)

% Attempt to Solve
% [x,fval,exitflag,info] = solve(Opt,x0)

clf;
opti_wikiplot('riemann');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2nlpa']);
end
clf;


% Ex6b

% Fitting Functions
ffun = @(x) [x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
             x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))];
% Objective
fun = @(x) norm(ffun(x));
% Bounds
lb = [-4;-4];
ub = [4;4];

% Starting Guess
x0 = [0;-3];

% Options
opts = optiset('solver','ipopt');

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)

% Attempt to Solve
[x,fval,exitflag,info] = solve(Opt,x0)
clf;
opti_wikiplot('wolfram');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2nlpb']);
end
clf;

%% GNLP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/GNLP

%Ex1
% Objective
fun = @(x) norm([x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
                 x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))]);

% Bounds
lb = [-4;-4]; 
ub = [4;4];

% Initial Solution
x0 = [-4;-4];

% Build optiprob structure (intermediate structure)
prob = optiprob('fun',fun,'bounds',lb,ub,'x0',x0);

% Setup solver options
opts1 = optiset('solver','nomad');
opts2 = optiset('solver','pswarm');
opts3 = optiset('solver','nlopt','solverOpts',nloptset('algorithm','GN_DIRECT'));
opts4 = optiset('solver','ipopt','warnings','off');

% Build OPTI objects
Opt1 = opti(prob,opts1); 
Opt2 = opti(prob,opts2); 
Opt3 = opti(prob,opts3); 
Opt4 = opti(prob,opts4);

% Solve the GNLP problems
[x1,fval1] = solve(Opt1,x0);
[x2,fval2] = solve(Opt2,x0);
[x3,fval3] = solve(Opt3,x0);
[x4,fval4] = solve(Opt4,x0);

%Ex2
% Objective
fun = @(x) x^6 - 2.08*x^5 + 0.4875*x^4 + 7.1*x^3 - 3.95*x^2 - x + 0.1;

% Bounds
lb = -1.5; 
ub = 1.5;

% Initial Solution
x0 = 0;

% Options
opts = optiset('solver','scip');

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'options',opts)

% Solve using SCIP
[x,fval,exitflag,info] = solve(Opt,x0)

% Plot Solution within Problem Bounds
plot(Opt,[lb ub])
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_gnlp_poly']);
end

%% MINLP
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/MINLP

%Ex1
% Objective
fun = @(x) -x(1) - x(2) - x(3);    

% Linear Constraints
A = [1 -1 0 0;
     1  0 1 1];
b = [0;2];

% Nonlinear Constraint
nlcon = @(x) (x(2) - 0.5)^2 + (x(3) - 0.5)^2; 
nlrhs = 0.25;
nle = -1; % -1 for <=, 0 for ==, +1 >=         

% Bounds
lb = [0;0;0;0];
ub = [1;Inf;Inf;5];

% Integer Constraints
xtype = 'BCCI';

%Initial Guess
x0 = [0;0;0;0];
% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'bounds',lb,ub,...
           'xtype',xtype)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)

% Ex2
% Objective
fun = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)))

% Bounds
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];

% Integer Constraints
xtype = 'IC';

%Initial Guess
x0 = [16;0];
% Options
opts = optiset('solver','nomad','display','iter')

% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'xtype',xtype,'options',opts)

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0)
clf;
opti_wikiplot('rastrigin');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'plot_ex2minlp']);
end
clf;

%% OPTIFIT
% https://www.inverseproblem.co.nz/OPTI/index.php/Probs/ModelFit

%Ex 1
% Fitting Data
xdata = (1:40)';
ydata = [5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237,...
         2.7002,2.8781, 2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791,...
         1.6686, 1.6232, 1.571, 1.6057,1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129,...
         1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,0.96518, 1.2129, 1.2003, 1.0743]';

% Automatically identify model structure and solve optimal parameters
ofit = optifit(xdata,ydata,'auto')

% Plot the Result
plot(ofit)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ofit_ex1']);
end

% Ex2
% Amount of Substrate
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; 
% Reaction Rate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; 

% Automatically identify model structure and solve optimal parameters
ofit = optifit(n,r,'auto')

% Plot the Result
plot(ofit)
h = legend;
set(h,'Location','northwest');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ofit_ex2']);
end

% Ex3
% Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure 
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate

% Model
Rxn = @(theta,xdata) theta(1)*xdata./(1+theta(2)*xdata);

% Fit Model Parameters
ofit = optifit(p,r,Rxn,[5e-3 2e-2])

% Plot the Result
plot(ofit)
h = legend;
set(h,'Location','northwest');
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ofit_ex3']);
end

%% Plot
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/Plots

%Ex1
% Objective (supplied with OPTI)
fun = @RiemannND;
% Bounds
lb = 0.5; ub = 2;
% Starting Guess
x0 = 1;

% Options
opts = optiset('solver','nomad');
% Create OPTI Object
Opt = opti('fun',fun,'bounds',lb,ub,'x0',x0,'options',opts)

% Plot Problem
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_1d1']);
end

% Plot Problem within Selected Bounds
plot(Opt,[lb ub])
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_1d2']);
end

% Plot Problem within Selected Bounds, and increased detail
plot(Opt,[lb ub],[],1000)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_1d3']);
end

% Solve the problem from improved initial guess
[x,fval] = solve(Opt,1.4)

% Plot the Problem again, this time include solution
plot(Opt,[lb,ub],[],1000)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_1d4']);
end

% 2D
% Objective
fun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2; 
% Constraints
lb = [-5;-5]; ub = [5;5];
% Initial Starting Guess
x0 = [0;0];

% Setup Options
opts = optiset('solver','lbfgsb');
% Create OPTI Object
Opt = opti('fun',fun,'x0',x0,'bounds',lb,ub,'options',opts)

% Plot
plot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_2d1']);
end

% Plot Problem +- 7 from x0
plot(Opt,7)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_2d2']);
end

% Plot Problem +- 7 from x0, and log(obj)
plot(Opt,7,true)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_2d3']);
end

% 3D
% QP Objective
H = eye(3);
f = -[2 3 1]';
% Linear Constraints
A = [1 1 1;3 -2 -3; 1 -3 2]; 
b = [1;1;1];
% Integer Constraints
xtype = 'CIC';

%Build & Solve
Opt = opti('qp',H,f,'ineq',A,b,'xtype',xtype)
try
    [x,fval,exitflag,info] = solve(Opt)
    % Plot Problem
    plot(Opt)
    set(gcf,'position',[510   344   688   558]); %BIGGER for higher D
    if(~isempty(figDir))
        print('-dpng','-r0',[figDir filesep 'ex_plot_3d1']);
    end
catch
end

% 5D
%Objective & Gradient
obj = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);
grad = @(x) [ x(1)*x(4) + x(4)*sum(x(1:3));
              x(1)*x(4);
              x(1)*x(4) + 1;
              x(1)*sum(x(1:3)) ];          
%Linear Constraints
lb = ones(4,1);
ub = 5*ones(4,1);
%Nonlinear Constraints
nlcon = @(x) [ prod(x);
               sum(x.^2)];
nljac = @(x) [ prod(x)./x';
                2*x' ];          
nlrhs = [25 40]';
nle = [1 0]'; % (>=, ==)
%Setup Options
opts = optiset('solver','ipopt','warnings','all','display','iter','derivCheck','on');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)
plot(Opt,2)
set(gcf,'position',[362   188   872   700]); %BIGGER for higher D
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_3d2']);
end
clf
set(gcf,'position',stdsize); %return to normal

% multi
% Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Linear Constraints
A = [-1 1];  b = -1;
Aeq = [1.1 1]; beq = 5; 
lb = [0;0]; ub = [4;4];
% Initial Guess
x0 = [2;2];

%Build & Solve using multiplot
Opt = opti('obj',obj,'bounds',lb,ub,'ineq',A,b,'eq',Aeq,beq)
[x,fval,ef,info] = multisolve(Opt,x0)

% Plot Problem using multiplot and take log(obj)
multiplot(Opt,1)
set(gcf,'position',[542   402   650   532]); %BIGGER for multi
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex_plot_m2d']);
end
set(gcf,'position',stdsize); %return to normal


%% Multi solve
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/MultiSolve

%Ex 1
% QP Objective
H = sparse([0 -1; -1 0]);
f = [0;0];
% Linear Constraints
lb = [-0.5;-0.5];
ub = [1;1];

% Options
opts = optiset('solver','clp','display','final');
% Build OPTI Object
Opt = opti('qp',H,f,'bounds',lb,ub,'options',opts)

% Solve Problem using MultiSolve
[x,fval,exitflag,info] = multisolve(Opt)

% Plot Problem + Search Area
multiplot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex1_multi1']);
end

% Ex2
% Solve Problem using MultiSolve and search 30^2 points
[x,fval,exitflag,info] = multisolve(Opt,[],30)

% Plot Problem + Search Area
multiplot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex1_multi2']);
end

% Ex3
% Solve Problem using MultiSolve [5x Phase 1, 5x Phase 2]
[x,fval,exitflag,info] = multisolve(Opt,[],[5 5])

% Plot Problem + Search Area
multiplot(Opt)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex1_multi3']);
end

% Ex4
% Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% Linear Constraints
A = [-1 1];  b = -1;
Aeq = [1.1 1]; beq = 5; 
lb = [0;0]; ub = [4;4];
% Initial Guess
x0 = [2;2];

% Options
opts = optiset('display','final');
% Build OPTI Object
Opt = opti('obj',obj,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'options',opts)

% Solve Problem using MultiSolve and Lowered Penalty
[x,fval,exitflag,info] = multisolve(Opt,x0,[],100)

% Plot Problem + Search Area
multiplot(Opt,1)
if(~isempty(figDir))
    print('-dpng','-r0',[figDir filesep 'ex4_multi1']);
end

%% SymBuilder
if (exist('syms.m','file'))
    % https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/SymBuilder
    % Create SymBuilder Object
    B = SymBuilder();
    % Create SymBuilder Object with suppressed command line output 
    Bs = SymBuilder(false);

    % Add as a String (scalar expression only)
    B.AddObj('-6*x1 -5*x2')

    % Add as a Symbolic Expression
    x = sym('x',[2 1]); %define symbolic variables
    Bs.AddObj(-[6;5]'*x)

    % Add as strings
    B.AddCon('x1 + 4*x2 <= 16');
    B.AddCon('6*x1 + 4*x2 <= 28')

    % Add as Symbolic Expressions with specified cl, cu
    Bs.AddCon(x(1) + 4*x(2),-Inf,16);
    Bs.AddCon(6*x(1) + 4*x(2),-Inf,28)

    % Method 1) Individually as Strings
    B.AddBound('0 <= x1 <= 10');
    B.AddBound('0 <= x2 <= 10')

    % Method 2) Vectorized, 'x' is recognised as belonging to all 'xn' variables
    B.AddBound('0 <= x <= 10');

    % Method 3) As a symbolic variable vector with numerical bounds
    Bs.AddBounds(x,[0;0],[10;10])

    % Build Optimization Problem
    Build(B)
    % Solve SymBuilder Optimization Problem
    [x,fval,ef,info] = Solve(B)

    % New SymBuilder Object
    B = SymBuilder();

    % Add Quadratic Objective
    B.AddObj('0.5*x1^2 + 0.5*x2^2 + 0.5*x3^2 - 2*x1 - 3*x2 - x3');

    % Add Linear Constraints
    B.AddCon('x1 + x2 + x3 <= 1');
    B.AddCon('3*x1 - 2*x2 - 3*x3 <= 1');
    B.AddCon('x1 - 3*x2 + 2*x2 <= 1');

    % Add Integer Constraint (also possible with Symbolic Variables)
    B.AddInteger('x2 = I');

    % Build Object
    Build(B)

    % Solve Problem
    try
        [x,fval,ef,info] = Solve(B)
    catch
    end

    % New SymBuilder Object
    B = SymBuilder();

    % Add Nonlinear Objective with 2 Constants
    B.AddObj('sin(pi*x1/a1)*cos(pi*x2/a2)');

    % Add Linear Constraints
    B.AddCon('-x1 + 2.5*x2 <= 1');
    B.AddCon('x1 + 2.5*x2 <= -15');

    % Add Integer Constraint (note also vectorized)
    B.AddInteger('x = I');

    % Declare Constants with Numerical Values
    B.AddConstant('a1',12);
    B.AddConstant('a2',16);

    % Build It
    Build(B)

    % Inspect resulting equations
    B.sobj
    % New SymBuilder Object
    B = SymBuilder();

    % Add Nonlinear Objective with 2 Intermediate Expressions
    B.AddObj('sin(e1)*cos(e2)');

    % Add Linear Constraints
    B.AddCon('-x1 + 2.5*x2 <= 1');
    B.AddCon('x1 + 2.5*x2 <= -15');

    % Add Integer Constraint
    B.AddInteger('x = I');

    % Add Intemediate Variable Expressions (with Constants for this example)
    B.AddExpression('e1 = pi*x1/a1');
    B.AddExpression('e2 = pi*x2/a2');

    % Define Constants
    B.AddConstant('a1',12);
    B.AddConstant('a2',16);

    % Add Now Build It
    Build(B)

    % New SymBuilder Object (verbose=false)
    B = SymBuilder(false);

    % Add Quadratic Objective
    B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');

    % Add Linear Constraints
    B.AddCon('x1 + 3*x2 = 5');
    B.AddCon('x3 + x4 - 2*x5 = 0');
    B.AddCon('x2 - x5 = 0');

    % Add Result Groups (the letter is arbitrary, any group name can be used)
    B.AddResultGroup('A','Fuel Usage [l]');
    B.AddResultGroup('B','Distance [km]');

    % Add Result Expressions (Group:Name, Variable or Expression)
    B.AddResultExp('A:Truck 1','x1');
    B.AddResultExp('A:Truck 2','x2');
    B.AddResultExp('B:Route 1','x3');
    B.AddResultExp('B:Route 2','x4');
    B.AddResultExp('B:Route 3','x5-x3'); %expressions with variables OK too

    % Build Object
    Build(B);

    % Solve
    Solve(B,[ 2.5 0.5 2 -1 0.5 ]);

    % Enter Problem
    B = SymBuilder();
    B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
    B.AddCon('x1^2 + 3*x2 = 4');
    B.AddCon('x3 + x4 - 2*x5 = 0');
    B.AddCon('x2 - x5 = 0');

    % Draft Model
    Draft(B)

    % Initial Guess
    x0 = [ 2.5 0.5 2 -1 0.5 ];

    % Generate C-Code Model and Solve Problem
    [x,fval,ef,info] = Solve(B,x0,symbset('cbmode','ccode'))

    % Generate C++ Code Model with CppAD and Solve Problem
    [x,fval,ef,info] = Solve(B,x0,symbset('cbmode','cppad'))

    % New SymBuilder Object
    B = SymBuilder(false);

    % Add Nonlinear Objective
    B.AddObj('sin(pi*x1/12)*cos(pi*x2/16)');

    % Add Linear Constraints
    B.AddCon('-x1 + 2.5*x2 <= 1');
    B.AddCon('x1 + 2.5*x2 <= -15');

    % Build Object
    Build(B);

    % Get Nonlinear Problem Data
    nlprob = GetNLProb(B)

    % Objective Function (min fun(x))
    fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);

    % Nonlinear Constraint Function (cl <= con(x) <= cu)
    con = @(x) [ prod(x);
                 sum(x.^2)];
    cl = [25;40];
    cu = [Inf;40];

    % Bounds + x0
    lb = ones(4,1);
    ub = 5*ones(4,1);
    x0 = [1 5 5 1]';

    % Build OPTI Object
    Opt = optisym(fun,x0,lb,ub,con,cl,cu)

    % Solve
    [x,fval,exitflag,info] = solve(Opt)

    B = SymBuilder();

    % Enter Objective
    B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');

    % Add Constraints
    B.AddCon('x1^2 + 3*x2 = 4');
    B.AddCon('x3 + x4 - 2*x5 = 0');
    B.AddCon('x2 - x5 = 0');

    Draft(B)
    Build(B)

    % Enter Problem
    B = SymBuilder();
    B.AddObj('(x1-x2)^2 + (x2+x3-2)^2 + (x4-1)^2 + (x5-1)^2');
    B.AddCon('x1^2 + 3*x2 = 4');
    B.AddCon('x3 + x4 - 2*x5 = 0');
    B.AddCon('x2 - x5 = 0');
    % Draft
    Draft(B)

    % Solve
    x0 = [ 2.5 0.5 2 -1 0.5 ];
    [x,fval,ef,info] = Solve(B,x0)

    % Solve without generating 2nd derivative callback
    [x,fval,ef,info] = Solve(B,x0,symbset('use2ndDerivs','no'))

    rehash
    clear symb_ccb
    clear symcb
    try
        delete('symb_ccb.mexw64')
    catch
    end
    try
        delete('symb_ccb.mexw32')
    catch    
    end
    try
       delete('symb_cb.m') 
    catch
    end
end

%% Two Cons
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/Cons

% Objective (f'x)
f = -[1 2 3]';

% Inequality Constraints (Ax <= b)
A = [-1  1  1;     
      1 -3  1];
b = [20;30];

% Equality Constraints (Aeqx = beq)
Aeq = [1 1 1];
beq = 40;

% Bounds (lb <= x <= ub)
lb = [0;0;0];
ub = [40;Inf;Inf];        

% Build OPTI Object
Opt = opti('f',f,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub)

% Solve Problem
[x,fval] = solve(Opt)

% Objective (f'x)
f = -[1 2 3]';

% Row Constraints (rl <= A*x <= ru)
A = [-1  1  1;     
      1 -3  1;
      1  1  1];
rl = [-Inf;-Inf;40]; 
ru = [20;30;40];

% Bounds (lb <= x <= ub)
lb = [0;0;0];
ub = [40;Inf;Inf];        

% Build OPTI Object
Opt = opti('f',f,'lin',A,rl,ru,'bounds',lb,ub)

% Solve Problem
[x,fval] = solve(Opt)

% Objective (fun(x))
fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);

% Mixed Nonlinear Constraints 
nlcon = @(x) [ prod(x); sum(x.^2) ];
nlrhs = [25;40];
nle = [1;0];      % [>=, =]

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

x0 = [1 5 5 1]';

% Build OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub)

% Solve Problem
[x,fval] = solve(Opt,x0)

% Objective (fun(x))
fun = @(x) x(1)*x(4)*sum(x(1:3)) + x(3);

% Row Nonlinear Constraints 
nlcon = @(x) [ prod(x); sum(x.^2) ];
cl = [25;40];
cu = [Inf;40]; 

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

x0 = [1 5 5 1]';

% Build OPTI Object
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub)

% Solve Problem
[x,fval] = solve(Opt,x0)

%% Low Level
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/LowLevel

% Ex1
% Load Large LP Test Problem
prob = coinRead('maros-r7.mps');

% Solve using OPTI Class
tic
solve(opti(prob,optiset('solver','clp')));
toc

% Solve using MATLAB Wrapper
tic
opti_clp([],prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub);
toc

% Solve using MEX Interface only
tic
clp([],prob.f,prob.A,prob.rl,prob.ru,prob.lb,prob.ub);
toc

% Ex2
% Fitting Function
fun = @(x) [100*(x(2)-x(1)^2); 1 - x(1)];

% Fitting Data
ydata = [0;0];

% Starting Guess
x0 = [-1.2;1];

% Call Low Level Interface
[x,fval] = opti_mkltrnls(fun,[],x0,ydata)

% Ex3
% Objective & Gradient
fun = @(x) log(1+x(1)^2) - x(2);
grad = @(x)[(2*x(1))/(x(1)^2+1), -1];

% Nonlinear Constraint, Jacobian & Structure
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2 - 4;
cl = 0; cu = 0;
nljac = @(x) sparse([4*x(1)*(x(1)^2+1),2*x(2)]);
jacstr = @() sparse([1 1]);        

% Starting Guess
x0 = [2;2];

% Build Function Structure
funcs.objective = fun;
funcs.gradient = grad;
funcs.constraints = nlcon;
funcs.jacobian = nljac;
funcs.jacobianstructure = jacstr;

% Build Options Structure
opts.lb = -1*Inf(2,1);
opts.ub = Inf(2,1);
opts.cl = cl;
opts.cu = cu;
opts.ipopt.hessian_approximation = 'limited-memory';

% Call IPOPT
[x,output] = ipopt(x0,funcs,opts)

% Build OPTI Object
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'x0',x0)

% Solve
[x,fval] = solve(Opt)


%% 1st Ders
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/Deriv1

% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);
% Initial Point
x0 = [1;2;3;4];

% Numerically Differentiate
df = mklJac(fun,x0)

% Automatic Differentiation
dfa = autoJac(fun,x0)

% Complex Step Differentiation
dfc = cstepJac(fun,x0)

% Symbolic Differentiation
if (exist('syms.m','file'))
    grad = symJac(fun)

    % Evaluate
    dfs = grad(x0)
end
    
% NLP Ex
% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Nonlinear Constraints 
nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
               x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
cl = [25;40];
cu = [Inf;40]; 

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

% Initial Guess
x0 = [1 5 5 1]';

% Gradient
grad = @(x) mklJac(fun,x);

% Jacobian
jac = @(x) mklJac(nlcon,x);

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,...
            'bounds',lb,ub,'x0',x0)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Gradient
grad = @(x) autoJac(fun,x);

% Jacobian
jac = @(x) autoJac(nlcon,x);

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,...
            'bounds',lb,ub,'x0',x0)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Gradient
grad = @(x) cstepJac(fun,x);

% Jacobian
jac = @(x) cstepJac(nlcon,x);

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,...
            'bounds',lb,ub,'x0',x0)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Gradient
grad = symJac(fun);

% Jacobian
jac = symJac(nlcon);

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,...
            'bounds',lb,ub,'x0',x0)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Gradient
grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3)), x(1)*x(4),...
             x(1)*x(4) + 1,  x(1)*sum(x(1:3))];

% Jacobian
jac = @(x) [prod(x')./x';
            2.*x'];

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,...
            'bounds',lb,ub,'x0',x0)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Sparsity
% Jacobian
jac = @(x) sparse([prod(x')./x';
                   2.*x']);
% Jacobian Structure
jacstr = @() sparse(ones(2,4));

% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Gradient
grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3)), x(1)*x(4),...
             x(1)*x(4) + 1,  x(1)*sum(x(1:3))];

% Nonlinear Constraints 
nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
               x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
cl = [25;40];
cu = [Inf;40]; 

% Jacobian
jac = @(x) sparse([prod(x')./x';
                   2.*x']);

% Jacobian Structure
jacstr = @() sparse(ones(2,4));

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

% Initial Guess
x0 = [1 5 5 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
            'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Options
opts = optiset('solver','ipopt','display','iter','derivCheck','on');

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
           'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)


%% 2nd Ders
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/Deriv2

% Objective
fun = @(x) x(1)*x(4)*(x(1) + x(2) + x(3)) + x(3);

% Gradient
grad = @(x) [x(1)*x(4) + x(4)*sum(x(1:3)), x(1)*x(4),...
             x(1)*x(4) + 1,  x(1)*sum(x(1:3))];

% Nonlinear Constraints 
nlcon = @(x) [ x(1)*x(2)*x(3)*x(4); 
               x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 ];
cl = [25;40];
cu = [Inf;40]; 

% Jacobian
jac = @(x) sparse([prod(x')./x';
                   2.*x']);

% Jacobian Structure
jacstr = @() sparse(ones(2,4));

% Hessian
H = @(x,sigma,lambda) sparse(sigma*[ 2*x(4)             0      0   0;
                                     x(4)               0      0   0;
                                     x(4)               0      0   0;
                                     2*x(1)+x(2)+x(3)  x(1)  x(1)  0 ] + ...
                     lambda(1)*[   0          0         0         0;
                                x(3)*x(4)     0         0         0;
                                x(2)*x(4) x(1)*x(4)     0         0;
                                x(2)*x(3) x(1)*x(3) x(1)*x(2)     0  ] + ...
                     lambda(2)*diag([2 2 2 2]));

% Hessian Structure
Hstr = @() sparse(tril(ones(4)));

% Bounds (lb <= x <= ub)
lb = ones(4,1);
ub = 5*ones(4,1);         

% Initial Guess
x0 = [1 5 5 1]';

% Options
opts = optiset('solver','ipopt','display','iter');

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
            'hess',H,'hstr',Hstr,'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)

% Options
opts = optiset('solver','ipopt','display','iter','derivCheck','on');

% Build OPTI Problem
Opt = opti('fun',fun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'jacstr',jacstr,...
            'hess',H,'hstr',Hstr,'bounds',lb,ub,'x0',x0,'options',opts)

% Solve NLP
[x,fval,exitflag,info] = solve(Opt)


%% MPS
% https://www.inverseproblem.co.nz/OPTI/index.php/File/MPS
 
prob = coinRead('testLP.mps')
 % Build an OPTI object of the returned problem 
OptMPS = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptMPS)

probq = coinRead('testQP2','qps',1)
Opt = opti(coinRead('testSOS2.mps'))

%% AMPL

prob = amplRead('diet.nl')
% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)

% Load an AMPL NLP
prob = amplRead('ch3.nl')

% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)

% Close ASL interface
asl('close')

 % Load an AMPL MILP
prob = amplRead('multmip1.nl')

% Build an OPTI object of the returned problem 
Opt = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(Opt)

try
    % Load an AMPL MIQCQP from Model File (but read as MINLP)
    prob = amplRead('trimlon.mod','trimlon2.dat',[],1)

    % Build an OPTI object of the returned problem 
    Opt = opti(prob)

    % Solve the resulting OPTI object
    [x,fval] = solve(Opt)

    % Load an AMPL MIQCQP from Model File
    prob = amplRead('trimlon.mod','trimlon2.dat')

    % Set SCIP as the solver
    opts = optiset('solver','scip','display','iter');

    % Build an OPTI object of the returned problem 
    Opt = opti(prob,opts)

    % Solve the resulting OPTI object
    [x,fval] = solve(Opt)
catch
end

% Load an AMPL NLP from Model File
prob = amplRead('hs100.nl')

% Examine constraint linearity
prob.conlin

%% SDPA
% https://www.inverseproblem.co.nz/OPTI/index.php/File/SDPA

prob = sdpRead('arch0.dat-s')
% Build an OPTI object of the returned problem 
OptSDPA = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptSDPA)

probd = sdpRead('sdpademo','SDPA',1)
prob = sdpRead('sdp_arch0.mat')
% Build an OPTI object of the returned problem 
OptSDMI = opti(prob)

% Solve the resulting OPTI object
[x,fval] = solve(OptSDMI)


%% WhiteBox Solvers
% https://www.inverseproblem.co.nz/OPTI/index.php/Advanced/WhiteBox
clc

% SCIP variable vector
x = scipvar(2,1);
% Test Function
fun = @(x) log(1 + x(1)^2) - x(2);

%Evaluate Function using SCIP variables
fun(x)

% Objective
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) [x(1) + x(2)^2;                    %cl <= nlcon(x) <= cu
              x(1)^2 + x(2);
              x(1)^2 + x(2)^2 - 1];
cl = [0;0;0];
cu = [Inf;Inf;Inf];

% Bounds
lb = [-0.5;-inf];                               %lb <= x <= ub
ub = [0.5;inf];

% Create OPTI Object
opts = optiset('solver','scip','display','iter');
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'opts',opts)

% Solve the NLP problem
x0 = [-2;1] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)

% Objective
fun = @(x) log(1 + x(1)^2) - x(2);    %Objective Function Vector (min f(x))

% Nonlinear Constraints
nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2; %Nonlinear Equality Constraint
cl = 4;
cu = 4;

% Create SCIP Options
sopts = scipset('scipopts',{'limits/solutions',1});

% Create OPTI Object
opts = optiset('solver','scip','display','iter','solverOpts',sopts);
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'ndec',2,'opts',opts)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)


% Objective
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

% Linear Constraints
A = [-1 1; 1 1];                                %rl <= A*x <= ru
rl = [-inf;5];
ru = [-1;5];

% Bounds
lb = [0;0];                                     %lb <= x <= ub
ub = [4;4];

% Integer Constraints
xtype = 'IC';

% Create OPTI Object
opts = optiset('solver','scip','display','iter');
Opt = opti('fun',fun,'lin',A,rl,ru,'bounds',lb,ub,'xtype',xtype,'opts',opts)

% Solve the NLP problem
x0 = [2;2] %Initial Guess
[x,fval,exitflag,info] = solve(Opt,x0)


% BARON variable vector
if (exist('barvec.m','file'))
    x = barvec(2,1);
    % Test Function
    fun = @(x) log(1 + x(1)^2) - x(2);

    %Evaluate Function using BARON variables
    fun(x)

    % Objective
    fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;    %Objective Function Vector (min f(x))

    % Nonlinear Constraints
    nlcon = @(x) [x(1) + x(2)^2;                    %cl <= nlcon(x) <= cu
                  x(1)^2 + x(2);
                  x(1)^2 + x(2)^2 - 1];
    cl = [0;0;0];
    cu = [Inf;Inf;Inf];

    % Bounds
    lb = [-0.5;-inf];                               %lb <= x <= ub
    ub = [0.5;inf];

    % Create OPTI Object
    opts = optiset('solver','baron','display','iter');
    Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'opts',opts)

    % Solve the NLP problem
    x0 = [-2;1] %Initial Guess
    [x,fval,exitflag,info] = solve(Opt,x0)


    % Objective
    fun = @(x) log(1 + x(1)^2) - x(2);    %Objective Function Vector (min f(x))

    % Nonlinear Constraints
    nlcon = @(x) (1 + x(1)^2)^2 + x(2)^2; %Nonlinear Equality Constraint
    cl = 4;
    cu = 4;

    % Create OPTI Object
    opts = optiset('solver','baron','display','iter');
    Opt = opti('fun',fun,'nl',nlcon,cl,cu,'ndec',2,'opts',opts)

    % Solve the NLP problem
    x0 = [2;2] %Initial Guess
    [x,fval,exitflag,info] = solve(Opt,x0)
end

%%
clc
fprintf('Done!\n');

