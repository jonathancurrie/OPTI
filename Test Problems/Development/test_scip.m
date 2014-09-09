%% Test Problems for SCIP Interface
% J.Currie 2012
clc
clear all

%% LP1 [-31.4]
clc
%Objective & Constraints
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]); 
rl = [-Inf;-Inf;-Inf];    
ru = [16;28;6];
lb = [];
ub = [];
opts = []; opts.display = 4; opts.maxiter = 1;
[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,[],[],[],[],opts)

%% LP2 [-97.5]
clc
%Objective & Constraints
f = -[1 2 3]';
A = sparse([-1,1,1; 1,-3,1; 1 1 1]);
rl = [-Inf;-Inf;40];
ru = [20,30,40]';
lb = [0;0;0];
ub = [40;inf;inf];

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub)

%% LP3 [2]
clc
%Objective & Constraints
f = [8,1]';
A = sparse([-1,-2;1,-4;3,-1;1,5;-1,1;-1,0;0,-1]); 
rl = -Inf(7,1);
ru = [-4,2,21,39,3,0,0]';
lb = [];
ub = [];

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub)

%% MILP1 [-29]
clc
%Objective & Constraints
f = -[6 5]';
A = sparse([1,4; 6,4; 2, -5]);
rl = -Inf(3,1);
ru = [16;28;6];  
lb = [0;0];
ub = [10;10];
xtype = 'II';
opts.display = 4;

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% MILP2 [-122.5]
clc
%Objective & Constraints
f = -[1 2 3 1]'; 
A = sparse([-1 1 1 10; 1 -3 1 0; 0 1 0 -3.5]); 
rl = [-Inf;-Inf;0];
ru = [20;30;0];  
lb = [0 0 0 2];
ub = [40 inf inf 3];
xtype = 'CCCI';

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype)

%% MILP3 Infeasible
clc
%Objective & Constraints
f = -[6 5]';
A = sparse([4,1; 1.5,1; -2, 1; -0.2, 1]); 
rl = -Inf(4,1);
ru = [5.3;4.5;-2.5; 1.5];  
lb = [0;0];
ub = [10;10];
xtype = 'II';

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype)

%% MILP4 [-3]
clc
%Objective & Constraints
f = -[-1, 2]';
A = sparse([2, 1;-4, 4]);
rl = [-Inf;-Inf];
ru = [5, 5]';
xtype = 'II';

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype)

%% BILP1 [-5]
clc
%Objective & Constraints
f = -[6 5]';
A = sparse([-3,5; 6,4; 3, -5; -6, -4]); 
rl = -Inf(4,1);
ru = [6;9;1;3];  
lb = []; ub = [];
xtype = 'BB';
opts.display = 5;

[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% BILP2 [-14]
clc
%Objective & Constraints
f = -[9 5 6 4]';
A = sparse([6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1]);
rl = -Inf(4,1);
ru = [9; 1; 0; 0];
xtype = 'BBBB';

[x,f,e,i] = scip([],f,A,rl,ru,[],[],xtype,[],[],[],opts)

%% LP SOS 1 [-230]
clc
f = [-1 -1 -3 -2 -2]';
A = sparse([-1 -1 1 1 0;
     1 0 1 -3 0]);
rl = [-Inf;-Inf];
ru = [30;30];
lb = -Inf(5,1);
ub = [40;1;inf;inf;1];

sos = '12';
sosind = {[1:2]' [3:5]'};
soswt = {[1:2]' [3:5]'};

s.type = sos;
s.index = sosind;
s.weight = soswt;
opts.display = 4;

[x,fval,e,i] = scip([],f,A,rl,ru,lb,ub,[],s,[],[],opts)

%% LP SOS 2 [-90]
clc
f = [-1 -1 -3 -2 -2]';
A = sparse([-1 -1 1 1 0;
             1 0 1 -3 0]);
rl = [-Inf;-Inf];
ru = [30;30];

lb = zeros(5,1);
ub = [40;1;inf;inf;1];

sos = '1';
sosind = [1:5]';
soswt = [1:5]';

s.type = sos;
s.index = sosind;
s.weight = soswt;
opts.display = 4;

[x,fval,e,i] = scip([],f,A,rl,ru,lb,ub,[],s,[],[],opts)

%% Large MILP (requires OPTI Toolbox)
clc
prob = coinRead('testMILP4.mps');
opts = optiset('solver','scip','display','iter');
Opt = opti(prob,opts)

[x,f,e,i] = solve(Opt);

%% QP1 [-2.833]
clc
%Objective & Constraints
H = speye(3);
f = -[2 3 1]';
A = sparse([1 1 1;3 -2 -3; 1 -3 2]); 
rl = -Inf(3,1);
ru = [1;1;1];
lb = []; ub = [];

xtype = 'CCC';
opts.display = 5;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% QP2 [-8.222]
clc
%Objective & Constraints
H = sparse([1 -1; -1 2]);
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1]);
rl = -Inf(3,1);
ru = [2; 2; 3];
lb = []; ub = [];

xtype = 'CC';
opts.display = 5;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% QP3 [-6.4138]
clc
%Objective & Constraints
H = sparse([1 -1; -1 2]);
f = -[2 6]';
A = sparse([1 1; -1 2; 2 1; 1 1.5]);
rl = [-Inf -Inf -Inf 2]';
ru = [2; 2; 3; 2];
lb = [0;0];
ub = [10;10];

xtype = 'CC';
opts.display = 5;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% Nonconvex QP [-1]
clc
%Objective & Constraints
H = sparse([0 -1; -1 0]);
f = [0;0];
%Constraints
lb = [-0.5;-0.5];
ub = [1;1];

xtype = 'CC';
opts.display = 5;

[x,f,e,i] = scip(H,f,[],[],[],lb,ub,xtype,[],[],[],opts)

%% MIQP1 [-11.5]
clc
%Objective & Constraints
H = sparse([1 -1; -1 2]);
f = -[2 6]';
A = sparse([1 1; -1 2]);
rl = [-Inf; -Inf];
ru = [3; 3.5];
xtype = 'IC';
lb = [0;0];
ub = [10;10];

opts.display = 5;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% MIQP2 [-2.75]
clc
%Objective & Constraints
H = speye(3);
f = -[2 3 1]';
A = sparse([1 1 1;3 -2 -3; 1 -3 2]); 
rl = -Inf(3,1);
ru = [1;1;1];
lb = []; ub = [];
xtype = 'CIC';

opts.display = 5;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],[],[],opts)

%% QCQP1 [-0.4242]
clc
%Objective & Constraints
H = sparse([33 6    0;
            6  22   11.5;
            0  11.5 11]);
f = [-1;-2;-3];
A = sparse([-1 1 1; 1 -3 1]);
rl = -Inf(2,1);
ru = [20;30];
Q = speye(3);
l = [0;0;0];
r = 1;
lb = [0;0;0];
ub = [40;inf;inf];

qc.Q = Q;
qc.l = l;
qc.qrl = -Inf(size(r));
qc.qru = r;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,[],[],qc,[],opts)

%% QCQP2 [-3.5]
clc
%Objective & Constraints
H = speye(2);
f = [-2;-2];
A = sparse([-1 1; 1 3]);
rl = -Inf(2,1);
ru = [2;5];
Q = sparse([1 0;0 1]);
l = [0;-2];
r = 1;
lb = [0;0];
ub = [40;inf];

qc.Q = Q;
qc.l = l;
qc.qrl = -Inf(size(r));
qc.qru = r;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,[],[],qc,[],opts)

%% QCQP3 [-1.7394]
clc
H = sparse(zeros(2));
f = [-2;-2];
A = sparse([-1 1; 1 3]);
rl = -Inf(2,1);
ru = [2;5];
Q = {sparse([1 0;0 1]);sparse([1 0;0 1])};
l = [0 2;2 -2];
r = [1 1];
lb = [0;0];
ub = [40;inf];

qc.Q = Q;
qc.l = l;
qc.qrl = -Inf(size(r));
qc.qru = r;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,[],[],qc,[],opts)

%% MIQCQP1 [-0.0942]
clc
%Objective & Constraints
H = sparse([33 6    0;
            6  22   11.5;
            0  11.5 11]);
f = [-1;-2;-3];
A = sparse([-1 1 1; 1 -3 1]);
rl = -Inf(2,1);
ru = [20;30];
Q = speye(3);
l = [0;0;0];
r = 1;
lb = [0;0;0];
ub = [40;inf;inf];
xtype = 'CCI';

qc.Q = Q;
qc.l = l;
qc.qrl = -Inf(size(r));
qc.qru = r;

[x,f,e,i] = scip(H,f,A,rl,ru,lb,ub,xtype,[],qc,[],opts)

%% Testing SCIP NONLINEAR EXPRESSIONS VIA TEST INTERFACE
clear all
clc

%% T1
clc
nl = @(x) 2*x;
x0 = -1;

opti_scipnltest(nl,x0)

%% T2
clc
nl = @(x) x*2;
x0 = -1;

opti_scipnltest(nl,x0)

%% T2
clc
nl = @(x) 2-2*x;
x0 = -1;
x = scipvar(1);
nl(x)

opti_scipnltest(nl,x0)

%% T3
clc
nl = @(x) 2+x*3;
x0 = -1;

opti_scipnltest(nl,x0)

%% T4
clc
nl = @(x) x*3+2;
x0 = -1;

opti_scipnltest(nl,x0)

%% T5
clc
nl = @(x) 2-x*3;
x0 = -1;

opti_scipnltest(nl,x0)

%% T5
clc
nl = @(x) x*3 - 2;
x0 = -1;

opti_scipnltest(nl,x0)

%% T6
clc
nl = @(x) 2/x*3;
x0 = -1;

opti_scipnltest(nl,x0)

%% T7
clc
nl = @(x) 2/(x*3);
x0 = -1;

opti_scipnltest(nl,x0)

%% T8
clc
nl = @(x) (2/x)*3;
x0 = -1;

opti_scipnltest(nl,x0)

%% T9
clc
nl = @(x) x*3 / (2*x);
x0 = -1;

opti_scipnltest(nl,x0)

%% T10
clc
nl = @(x) (1-x) / (1+x)*2;
x0 = -2;

opti_scipnltest(nl,x0)

%% T11
clc
nl = @(x) (1-x) / ((1+x)*2);
x0 = -3;

opti_scipnltest(nl,x0)

%% T12
clc
nl = @(x) 1-(x+2);
x0 = -1;

opti_scipnltest(nl,x0)

%% T
clc
nl = @(x) 3*x/4;
x0 = -1;

opti_scipnltest(nl,x0)

%% T13
clc
nl = @(x) 2*x(1) - 3*x(2);
x0 = [-1;-1];

opti_scipnltest(nl,x0)

%% T13
clear all
clc
nl = @(x) 3*x(1)*x(2) + 2*x(1) - 3*x(2)/4;
x0 = [-1;-1];
x = scipvar(2,1);
nl(x)

opti_scipnltest(nl,x0)

%% T14
clc
nl = @(x) 0.5*x(2) - x(1)*2;
x0 = [-1;-1];

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) [3*x(1)*x(2) + 2*x(1) - 3*x(2)/4
           2*x(1)*x(2)];
x0 = [-1;-1];

opti_scipnltest(nl,x0)


%%
clc
nl = @(x) [x(1)-x(2)^2;
           x(1)*x(2) + 2*x(1) - 3*x(2)/4;
           10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1)];
x0 = [-3  -1  -3  -1]';

x = scipvar(4,1);
nl
n = nl(x);
n(1)
n(2)
n(3)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) [10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1)];
x0 = [-0.5  -2  -3  -1]';

x = scipvar(4,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);
x0 = [-3  -1  -3  -1]';

x = scipvar(4,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
x0 = [-1;-2;-3;-4;-5];           

x = scipvar(5,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) -x(1) - x(2) - x(3);
x0 = [-1;-2;-3];

x = scipvar(3,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) x(2) + (x(2)-x(1))^2;
x0 = [-1;-2];

x = scipvar(2,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)
%%
clc
nl = @(x) (x(1)*x(2))^((x(2)-x(1))^2);
x0 = randn(2,1);

x = scipvar(2,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) x(3) - 2*(x(4) + x(5)) + 3;
x0 = [-1;-2;-3;-4;-5];

x = scipvar(5,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
nl = @(x) exp(x(1));

x0 = [0;1.05];
x = scipvar(2,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
fun = @(x) -x(1);
nl = @(x) [x(1) - exp(x(2));
              x(3) - exp(x(2))];

x0 = [0;1.05;2.9];
x = scipvar(5,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
clear all
nl = @(x)  x(1) - (5.1*x(4)*(1.2*x(1)/4.5 - 3.4));

x0 = randn(27,1);%(0:26)';
x = scipvar(27,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
clear all
nl = @(x) x(19)*((7786463903543161*x(2))/9007199254740992 + 31189022695929/137438953472) + x(3)*x(24) + x(6)*((6881844292114419*x(1))/9007199254740992 + 6863778725243905/17592186044416) - x(3)*(x(25) + x(26) + x(27)) + x(4)*(x(1) - ((log((5103105382857231*x(4)*((1518110687590409*x(1))/4503599627370496 - 1225674772364983/2199023255552))/147573952589676412928 + 1/10) + 5886760970214945/1125899906842624)*((1518110687590409*x(1))/4503599627370496 - 1225674772364983/2199023255552))/7) + x(15)*((115734561392307911515507040800193*x(1))/162259276829213363391578010288128 + 10201874622980696614228209270927/19807040628566084398385987584)

x0 = randn(27,1);%(0:26)';
x = scipvar(27,1);
nl
n = nl(x);
n(1)

opti_scipnltest(nl,x0)

%%
clc
% clear all
nl = @(x) x(1);
x0 = randn(1,1);
opti_scipnltest(nl,x0)


%% ERROR ON WIN32 (NO IDEA WHY)
clc
clear all
fun = @(x) -x(1)*x(2)*x(3);
nlcon = @(x) -x(1)^2 - 2*x(2)^2 - 4*x(3)^2 + 48;
cl = 0;
cu = inf;
x0 = [1;1;1];
fmin = -16*sqrt(2);
opts = optiset('display','iter');

lb = [0;0;0];
ub = [5;5;5];
% opti_scipnltest(nlcon,x0)

x = scipvar(3,1);
nlcon(x)

[x,fval,ef] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts)


%% UNO Rosenbrock
clc
%Objective
obj = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
%Setup Options
opts = optiset('solver','nlopt');
%Build & Solve
Opt = opti('obj',obj,'ndec',2,'options',opts)
x0 = [0 0]';
[x,fval,exitflag,info]= solve(Opt,x0)

%%
[x,fval,exitflag,info] = opti_scipnl(obj,[],[],[],[],[],[],[],[],[],x0)

%% NLP1 Hock & Schittkowski #71
clc
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
opts = optiset('solver','ipopt','warnings','on','display','iter');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'bounds',lb,ub,'options',opts)
x0 = [1 5 5 1]';
[x,fval,exitflag,info]= solve(Opt,x0)

info.Lambda

%%
cl = [25; 40];
cu = [Inf; 40];

[x,fval,exitflag,info] = opti_scipnl(obj,[],[],[],lb,ub,nlcon,cl,cu,[],x0)

%% NLP2 Hock & Schittkowski #38
clc
%Objective & Gradient
obj = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + 90*(x(4)-x(3)^2)^2 + (1-x(3))^2 + ...
      10.1*(x(2)-1)^2 + 10.1*(x(4)-1)^2 + 19.8*(x(2)-1)*(x(4)-1);  
grad = @(x) [ -400*x(1)*(x(2)-x(1)^2) - 2*(1-x(1));
              200*(x(2)-x(1)^2) + 20.2*(x(2)-1) + 19.8*(x(4)-1);
              -360*x(3)*(x(4)-x(3)^2) - 2*(1-x(3));
              180*(x(4)-x(3)^2) + 20.2*(x(4)-1) + 19.8*(x(2)-1)];
%Hessian (must not be tril for matlab!) & Hessian Structure       
hess = @(x) sparse(  [ 1200*x(1)^2-400*x(2)+2  -400*x(1)       0                          0
                       -400*x(1)               220.2           0                          19.8
                        0                      0               1080*x(3)^2- 360*x(4) + 2  -360*x(3)
                        0                      19.8            -360*x(3)                   200.2 ]);  
Hstr = @() sparse([ 1  0  0  0 
                    1  1  0  0
                    0  0  1  0
                    0  1  1  1 ]);                               
%Constraints
lb = [-10 -10 -10 -10]';
ub = [10 10 10 10]';
%Setup Options
opts = optiset('solver','ipopt');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'hess',hess,'Hstr',Hstr,'bounds',lb,ub,'options',opts)
x0 = [-3  -1  -3  -1]';
[x,fval,exitflag,info]= solve(Opt,x0)

%%
[x,fval,exitflag,info] = opti_scipnl(obj,[],[],[],lb,ub,[],[],[],[],x0,opts)

%% NLP3 Hock & Schittkowski #51
clc
%Objective & Gradient
obj = @(x) (x(1) - x(2))^2 + (x(2) + x(3) - 2)^2 + (x(4) - 1)^2 + (x(5) - 1)^2;
grad = @(x) 2*[ x(1) - x(2);
                x(2) + x(3) - 2 - x(1) + x(2);
                x(2) + x(3) - 2;
                x(4) - 1;
                x(5) - 1 ];
%Nonlinear Constraints & Jacobian Structure
nlcon = @(x) [ x(1) + 3*x(2);
               x(3) + x(4) - 2*x(5);
               x(2) - x(5) ];
nljac = @(x) sparse([ 1  3  0  0  0;
	                  0  0  1  1 -2;
	                  0  1  0  0 -1 ]);
nljacstr = @() sparse([1 1 0 0 0;
                       0 0 1 1 1;
                       0 1 0 0 1]);
nlrhs = [4 0 0]';
nle = [0 0 0]';
%Setup Options
opts = optiset('solver','nlopt');
%Build & Solve
Opt = opti('obj',obj,'grad',grad,'nlmix',nlcon,nlrhs,nle,'nljac',nljac,'nljacstr',nljacstr,'options',opts)
x0 = [ 2.5 0.5 2 -1 0.5 ];
[x,fval,exitflag,info]= solve(Opt,x0)

%%
cl = [4;0;0]; cu = cl;
[x,fval,exitflag,info] = opti_scipnl(obj,[],[],[],[],[],nlcon,cl,cu,[],x0,opts)

%% NLP4 Hock & Schittkowski #104 (not working)
clc
%Objective
obj = @(x) 0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
%Constraints
nlcon = @(x) [ 1 - 0.588*x(5)*x(7) - 0.1*x(1);
               1 - 0.588*x(6)*x(8) - 0.1*x(1) - 0.1*x(2);
               1 - 4*x(3)*x(5)^(-1) - 2*x(3)^(-0.71)*x(5)^(-1) - 0.588*x(3)^(-1.3)*x(7);
               1 - 4*x(4)*x(6)^(-1) - 2*x(4)^(-0.71)*x(6)^(-1) - 0.588*x(4)^(-1.3)*x(8);
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2);
               0.4*x(1)^0.67 * x(7)^(-0.67) + 0.4*x(2)^0.67 * x(8)^(-0.67) + 10 - x(1) - x(2)];
nlrhs = [ 0 0 0 0 1 4.2]';
nle = [1 1 1 1 1 -1]';
lb = 0.1*ones(8,1);
ub = 10*ones(8,1);
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
x0 = [6,3,.4,.2,6,6,1,.5];
Opt = opti('obj',obj,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info]= solve(Opt,x0)

%%
cl = [0;0;0;0;1;-Inf];
cu = [Inf;Inf;Inf;Inf;Inf;4.2];
[x,fval,exitflag,info] = opti_scipnl(obj,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts)

%% Can't solve...

fun = @(x) (x(1)-x(2))^2 + (x(2)-x(3))^3 + (x(3)-x(4))^4 + (x(4)-x(5))^4;
nlcon = @(x) [x(1) + x(2)^2 + x(3)^3 - 3;
              x(2) - x(3)^2 + x(4) - 1;
              x(1)*x(5) - 1];
cl = [0;0;0];
cu = [0;0;0];
lb = -10*ones(5,1); ub = -1*lb;
x0 = [2;sqrt(2);-1;2*-sqrt(2);0.5];
% fmin(i) = 0;
opts = optiset('display','iter','maxnodes',1e5);

[x,fval,ef,info] = opti_scipnl(fun,[],[],[],lb,ub,nlcon,cl,cu,[],x0,opts)

%% NLP testing
clc
clear
fun = @(x) -2*x(1) + x(2);
f = [-2 1];

nlcon = @(x) 2.5*x(1)*x(2) - 0.2*x(2) + 0.5*x(1) + x(1)*x(1)*x(1);
nlrhs = 1;
nle = 0;

SCIP_NUM = 0;
SCIP_VAR = 1;
SCIP_MUL = 3;
SCIP_DIV = 4;
SCIP_ADD = 5;
SCIP_SUB = 6;

nl.instr = [SCIP_NUM; 2.5;
            SCIP_VAR; 0;
            SCIP_MUL; NaN
            SCIP_VAR; 1;
            SCIP_MUL; NaN;
            SCIP_NUM; 0.2;
            SCIP_VAR; 1;
            SCIP_MUL; NaN;
            SCIP_SUB; NaN;
            SCIP_NUM; 0.5;
            SCIP_VAR; 0;
            SCIP_MUL; NaN;
            SCIP_ADD; NaN;
            SCIP_VAR; 0;
            SCIP_VAR; 0;
            SCIP_MUL; NaN;
            SCIP_VAR; 0;
            SCIP_MUL; NaN;
            SCIP_ADD; NaN];
nl.cl = 1;
nl.cu = 1;

A = []; rl = []; ru = [];
lb = [-1;-1];
ub = [1;1];
xtype = 'CC';
sos =[];
opts.display = 5;
      
[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype,sos,[],nl,opts)
     
%

Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'x0',[1;1],'options',optiset('solver','ipopt'));

[x,f,e,i] = solve(Opt)
plot(Opt)

%% baron

[x,f,e,i] = baron(fun,[],[],[],lb,ub,nlcon,1,1)



%% NLP testing
clc
clear
fun = @(x) -2*x(1) + x(2);
f = [-2 1];

nlcon = @(x) 0.5*x(2) - x(1)*2;
% nlcon = @(x) x(2) - (2.5*x(1));
nlrhs = 1;
nle = 0;

SCIP_NUM = 0;
SCIP_VAR = 1;
SCIP_MUL = 3;
SCIP_DIV = 4;
SCIP_ADD = 5;
SCIP_SUB = 6;

nl.instr = [SCIP_VAR; 0;
            SCIP_NUM; 2;
            SCIP_MUL; NaN;
            SCIP_NUM; 0.5;
            SCIP_VAR; 1;
            SCIP_MUL; NaN;
            SCIP_SUB; 1];
nl.cl = 1;
nl.cu = 1;

A = []; rl = []; ru = [];
lb = [-1;-1];
ub = [1;1];
xtype = 'CC';
x0 = [1;1];
sos =[];
opts.display = 5;
opts.x0 = x0;
opts.nlcon_val = nlcon(x0);
      
[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype,sos,[],nl,opts)
      
%

Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'x0',x0,'options',optiset('solver','matlab'));

[x,f,e,i] = solve(Opt)
plot(Opt)

%% baron

[x,f,e,i] = baron(fun,[],[],[],lb,ub,nlcon,1,1)

%% NLP testing
clc
clear
fun = @(x) -2*x(1) + x(2);
f = [-2 1];

nlcon = @(x) 0.3*x(1).^(-2.7);
% nlcon = @(x) x(2) - (2.5*x(1));
nlrhs = 1;
nle = 0;

SCIP_NUM = 0;
SCIP_VAR = 1;
SCIP_MUL = 3;
SCIP_DIV = 4;
SCIP_ADD = 5;
SCIP_SUB = 6;
SCIP_SQU = 7;
SCIP_SQR = 8;
SCIP_POW = 9;
SCIP_EXP = 10;
SCIP_LOG = 11;

nl.instr = [SCIP_VAR; 0;
            SCIP_NUM; -2.7;
            SCIP_POW; NaN;
            SCIP_NUM; 0.3;
            SCIP_MUL; NaN];
nl.cl = 1;
nl.cu = 1;

A = []; rl = []; ru = [];
lb = [-1;-1];
ub = [1;1];
xtype = 'CC';
sos =[];
opts.display = 5;
      
[x,f,e,i] = scip([],f,A,rl,ru,lb,ub,xtype,sos,[],nl,opts)
      
%

Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'x0',[1;1],'options',optiset('solver','matlab'));

[x,f,e,i] = solve(Opt)
plot(Opt)

%% baron

[x,f,e,i] = baron(fun,[],[],[],lb,ub,nlcon,1,1)

   