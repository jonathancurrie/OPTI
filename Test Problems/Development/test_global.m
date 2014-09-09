%% A Collection of Global Optimization Problems
% N. Sahinidis (CMU)
%#ok<*NASGU,*ASGLU,*NOPTS>

%% St_e01 [x = 6,0.6667, fval = -6.6667]
clc
%Objective
fun = @(x) -x(1) - x(2);
%Nonlinear Constraints
nlcon = @(x) x(1)*x(2);
cl = -Inf;
cu = 4;
%Bounds
lb = [0;0];
ub = [6;4];
x0 = [1,1]';

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Plot
plot(Opt)
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e02 [x = 6.2934, 3.8218, 201.590, fval = 201.1591]
clc
%Objective
fun = @(x) x(3);
%Nonlinear Constraints
nlcon = @(x) [ 30*x(1) - 6*x(1)*x(1) - x(3);
               20*x(2) - 12*x(2)*x(2) - x(3);
               0.5*(x(1) + x(2))^2 - x(3) ];
cl = [-250;-300;-150];
cu = [-250;-300;-150];
%Bounds
lb = [0;0;0];
ub = [9.422;5.9023;267.417085245];
x0 = [1,1,1];

%Setup Options
opts = optiset('solver','gmatlab','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e03 [fval = -1.1613e3]
clc
%Objective
fun = @(x) -0.063*x(4)*x(7) + 5.04*x(1) + 0.035*x(2) + 10*x(3) + 3.36*x(5);
%Linear Constraints
A = [1 0 0 -1.22 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0.222;
     0 0 0 0 0 0 3 0 0 -1];
rl = [0;35.82;133];
ru = [0;35.82;133];                        
%Nonlinear Constraints
nlcon = @(x) [ 0.038*x(8)^2 - 1.098*x(8) - 0.325*x(6) + x(7);
               x(4)*x(9)*x(6) + 1000*x(3)*x(6) - 98000*x(3);
               -x(1)*x(8) + x(2) + x(5);
               0.13167*x(8)*x(1) + 1.12*x(1) - 0.00667*x(8)^2*x(1) - x(4)];
cl = [57.425;0;0;0];
cu = [57.425;0;0;Inf];
%Bounds
lb = [1;1;0;1;0;85;90;3;1.2;145];
ub = [2000;16000;120;5000;2000;93;95;12;4;162];
%Starting point
x0 = [1;1;1;1;1;85;90;3;1.2;145];

%Setup Options
opts = optiset('solver','gmatlab','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'lin',A,rl,ru,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e04 [x = 0, 94.1779, 0.0001, 80, fval = 5.1949e3]
clc
%Objective
fun = @(x) 400*x(1)^0.9 + 22*(-14.7 + x(2))^1.2 + x(3) + 1000;
%Nonlinear Constraints
nlcon = @(x) [ x(3)*x(1) + 144*x(4);
               -exp(11.86 - 3950/(460 + x(4))) + x(2) ];
cl = [11520;0];
cu = [Inf;0];
%Bounds
lb = [0;14.7;0;-459.67];
ub = [15.1;94.2;5371;80];
%Starting point
x0 = [1;14.7;1;1];

%Setup Options
opts = optiset('solver','gmatlab','display','iter','solverOpts',psoptimset('MaxFunEvals',10e4)); 
%Build & Solve
Opt = opti('fun',fun,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e05 [fval = 7.0492e3]
clc
%Objective
fun = @(x) x(1) + x(2) + x(3);
%Linear Constraints
A = [0 0 -4e3 0 -1e5];
rl = -5e7;
ru = -5e7;
%Nonlinear Constraints
nlcon = @(x) [ 1e5*x(4) - 120*x(1)*(300-x(4));
               1e5*x(5) - 80*x(2)*(400-x(5)) - 1e5*x(4)];
cl = [1e7;0];
cu = [1e7;0];
%Bounds
lb = [0;0;0;100;100];
ub = [15834;36250;10000;300;400];
%Starting point
x0 = [1;1;1;100;100];

%Setup Options
opts = optiset('solver','gmatlab','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'lin',A,rl,ru,'nl',nlcon,cl,cu,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e08 [x = 0.1294, 0.4830, fval = 0.7418]
clc
%Objective
fun = @(x) 2*x(1)+x(2);
%Constraints
nlcon = @(x) [-16*x(1)*x(2); (-4*x(1)^2) - 4*x(2)^2];
nlrhs = [-1;-1];
nle = [-1;-1];
lb = [0;0];
ub = [1;1];
x0 = [0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

% %% St_e09 [x = 0.5, 0.5, fval = -0.5]
% clc
% %Objective
% fun = @(x) -2*x(1)*x(2);
% %Constraints
% nlcon = @(x) 4*x(1)*x(2) + 2*x(1) + 2*x(2);
% nlrhs = 3;
% nle = -1;
% lb = [0;0];
% ub = [1;1];
% x0 = [0;0];
% 
% %Setup Options
% opts = optiset('solver','nomad','display','iter'); 
% %Build & Solve
% Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
% [x,fval,exitflag,info] = solve(Opt,x0)  
% %Check Solution
% [ok,msg] = checkSol(Opt)

%% St_e17 [x = 8.1651, 7.5685, fval = 376.2889]
clc
%Objective
fun = @(x) 29.4*x(1) + 18*x(2);
%Constraints
nlcon = @(x) -(x(1) - 0.2458*x(1)^2/x(2));
nlrhs = -6;
nle = -1;
lb = [0;1e-5];
ub = [115.8;30];
x0 = [0;1e-5];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e18 [x = -1.4142, -1.4142, fval = -2.8284]
clc
%Objective
fun = @(x) x(1) + x(2);
%Constraints
nlcon = @(x) [-(x(1)^2 + x(2)^2); x(1)^2 + x(2)^2; -x(1) + x(2); x(1) - x(2)];
nlrhs = [-1;4;1;1];
nle = [-1;-1;-1;-1];
lb = [-2;-2];
ub = [2;2];
x0 = [0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% St_e19 [x = -3.1736, 1.7245, fval = -118.7049]
clc
%Objective
fun = @(x) x(1)^4 - 14*x(1)^2 + 24*x(1) - x(2)^2;
%Constraints
nlcon = @(x) [-x(1) + x(2); (-x(1)^2) - 2*x(1) + x(2)];
nlrhs = [8;-2];
nle = [-1;-1];
lb = [-8;0];
ub = [10;10];
x0 = [0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter','maxtime',5); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% MINLP1 
clc
%Objective
fun = @(x) (x(1) - 5)^2 + x(2)^2 - 25;
%Constraints
nlcon = @(x) x(1)^2 - x(2) + 0.5;
nlrhs = 0; nle = -1;    
xtype = 'II';
x0 = [0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% MINLP nvs03
clc
%Objective
fun = @(x) (-8 + x(1))^2 + (-2 + x(2))^2;
%Constraints
nlcon = @(x) [-0.1*x(1)^2 + x(2); -1/3*x(1) - x(2)];
nlrhs = [0;-4.5]; nle = [1;1];    
lb = [0;0]; ub = [200;200];
xtype = 'II';
x0 = [0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% MINLP nvs04
clc
%Objective
fun = @(x) 100*(0.5 - (0.6 + x(1))^2 + x(2))^2 + (0.4 - x(1))^2;
%Constraints  
lb = [0;0]; ub = [200;200];
xtype = 'II';
x0 = [4;4];

%Setup Options
opts = optiset('solver','nomad','display','iter','solverOpts',nomadset('vns_search',0.9)); 
%Build & Solve
Opt = opti('fun',fun,'bounds',lb,ub,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)
plot(Opt,[],1);

%% MINLP nvs06
clc
%Objective
fun = @(x) 0.1*(x(1)^2 + (1 + x(2)^2)/x(1)^2 + (100 + x(1)^2*x(2)^2)/(x(1)*x(2))^4) + 1.2;
%Constraints  
lb = [1;1]; ub = [200;200];
xtype = 'II';
x0 = [1;2];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'bounds',lb,ub,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)
plot(Opt,2,1);

%% MINLP nvs07
clc
%Objective
fun = @(x) 2*x(2)^2 + x(1) + 5*x(3);
%Constraints
A = [1 0 -1]; b = 2.66; e = 1;
nlcon = @(x) x(3)^2*x(2) + 5*x(3) + 3*x(1);
nlrhs = 10; nle = 1;    
lb = [0;0;0]; ub = [200;200;200];
xtype = 'III';
x0 = [0;0;0];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'mix',A,b,e,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% MINLP nvs08
clc
%Objective
fun = @(x) (-3 + x(1))^2 + (-2 + x(2))^2 + (4 + x(3))^2;
%Constraints
nlcon = @(x) [x(3)^0.5 + x(1) + 2*x(2); 
              0.240038406144983*x(1)^2 - x(2) + 0.255036980362153*x(3); 
              x(2)^2 - 1/(x(3)^3*x(3)^0.5) - 4*x(1)];
nlrhs = [10;-3;-12]; nle = [1;1;1];    
lb = [0;0;0.001]; ub = [200;200;200];
xtype = 'IIC';
x0 = [0;0;0.001];

%Setup Options
opts = optiset('solver','nomad','display','iter'); 
%Build & Solve
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'bounds',lb,ub,'ctype',xtype,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)  
%Check Solution
[ok,msg] = checkSol(Opt)

%% MINLP MATLAB GA Example
clc
%Objective
obj = @(x) 20 + x(1)^2 + x(2)^2 - 10*(cos(2*pi*x(1)) + cos(2*pi*x(2)));
%Constraints
lb = [5*pi;-20*pi];
ub = [20*pi;-4*pi];
%Setup Options
opts = optiset('solver','nomad','display','iter');
% Solve
Opt = opti('obj',obj,'bounds',lb,ub,'int','IC','options',opts)
x0 = [60;0];
[x,fval,ef,info] = solve(Opt,x0)