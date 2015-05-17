%% Testing SCIP Option Setting Methods
% 17/5/15
clc
clear
% Default Problem
a = 100;
fun = @(x) a*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
lb = [-1.5; -1.5]; ub = [2;2]; x0 = [-2; 1];    
%Convert to Ins List
zs = zeros(size(x0));
x = scipvar(size(x0));
fs = fun(x); nl.obj_instr = fs.ins; nl.obj_val = fun(x0); nl.x0 = x0;
opts = optiset('solver','scip');

%% No scipopts Field
sopts = [];
sopts.maxiter = 1;
sopts.display = 3;
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Not a cell array
sopts.scipopts = 'a';
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Cell array of one column
sopts.scipopts = {'a';1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Empty Name
sopts.scipopts = {'',1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Empty Val (will skip it)
sopts.scipopts = {'a',[]};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Unknown Option
sopts.scipopts = {'limits/times',1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Wrong Datatype for Option
sopts.scipopts = {'limits/time','a'};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Outside Option Limits
sopts.scipopts = {'limits/time',-1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Bool Option logical
sopts.scipopts = {'display/lpinfo',false};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Bool Option double
sopts.scipopts = {'display/lpinfo',1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Char Option not char
clc
sopts.scipopts = {'lp/initalgorithm',1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Char Option not valid char
clc
sopts.scipopts = {'lp/initalgorithm','a'};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% Char Option 
clc
sopts.scipopts = {'lp/initalgorithm','p'};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% String Option not string
clc
sopts.scipopts = {'nlp/solver',1};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% String Option unknown
clc
sopts.scipopts = {'nlp/solver','a'};
scip([],zs,[],[],[],lb,ub,[],[],[],nl,sopts)

%% opti_scipnl test
clc
fun = @(x) -x(1)*x(2)*x(3);
nlcon = @(x) -x(1)^2 - 2*x(2)^2 - 4*x(3)^2 + 48;
cl = 0;
cu = inf;
x0 = [1;1;1];
opts = optiset('display','iter','solverOpts',scipset('scipopts',{'nlp/disable',true}));
opti_scipnl(fun,[],[],[],[],[],nlcon,cl,cu,[],x0,opts);


%% opti test
clc
fun = @(x) (x(1)-x(2))^2 + (x(3)-1)^2 + (x(4)-1)^4 + (x(5)-1)^6;
nlcon = @(x) [x(1) + x(2) + x(3) + 4*x(4) - 7;
              x(3) + 5*x(5) - 6];
cl = [0;0];
cu = [0;0];
x0 = [10;7;2;-3;0.8];

sopts = scipset('scipopts',{'limits/solutions',1}); %find first feasible solution
opts = optiset('solver','scip','display','iter','solverOpts',sopts);
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'x0',x0,'opts',opts)

[x,f,e,i] = solve(Opt)





