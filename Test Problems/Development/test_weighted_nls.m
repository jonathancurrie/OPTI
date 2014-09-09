%% NLS1 Original
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Setup Options
opts = optiset('solver','matlab','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
confidence(Opt)
plot(Opt)

%% NLS1 Weighted [no grad]
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',[],'data',xdata,ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
[cInt,stats] = confidence(Opt)
plot(Opt)

%% NLS1 Weighted no xdata [no grad]
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Function
fun = @(x) x(1)*exp(x(2)*xdata);
grad = @(x) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',[],'ydata',ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted [with grad]
clc
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted no xdata [with grad]
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Function
fun = @(x) x(1)*exp(x(2)*xdata);
grad = @(x) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','lmder','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'ydata',ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP [no grad], MANUAL WEIGHTS
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
%Fitting Weights
weights = ones(size(xdata)); weights(end) = 1e3; 
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]'.*weights;
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata).*weights;
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',[],'data',xdata,ydata,'weights',[],'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP [no grad], AUTO WEIGHTS
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]';
%Fitting Weights
weights = ones(size(xdata)); weights(end) = 1e3; 
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',[],'data',xdata,ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP no xdata [no grad]
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Function
fun = @(x) x(1)*exp(x(2)*xdata);
grad = @(x) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',[],'ydata',ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP [WITH grad], MANUAL WEIGHTS
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
%Fitting Weights
weights = ones(size(xdata)); weights(end) = 1e3;
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]'.*weights;
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata).*weights;
grad = @(x,xdata) [ exp(x(2).*xdata).*weights, x(1).*xdata.*exp(x(2).*xdata).*weights];
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'weights',[],'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP [WITH grad], AUTO WEIGHTS
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5]';
%Fitting Weights
weights = ones(size(xdata)); weights(end) = 1e3;
%Function
fun = @(x,xdata) x(1)*exp(x(2)*xdata);
grad = @(x,xdata) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'data',xdata,ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS1 Weighted as NLP no xdata [WITH grad]
clc
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
%Function
fun = @(x) x(1)*exp(x(2)*xdata);
grad = @(x) [ exp(x(2).*xdata), x(1).*xdata.*exp(x(2).*xdata)];
%Fitting Weights
weights = ones(size(ydata)); weights(end) = 1e3;
%Setup Options
opts = optiset('solver','ipopt','display','iter');
%Build & Solve
Opt = opti('fun',fun,'grad',grad,'ydata',ydata,'weights',weights,'ndec',2,'options',opts)
x0 = [100; -1]; % Starting guess
[x,fval,exitflag,info] = solve(Opt,x0)
%Plot
plot(Opt)

%% NLS2
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Setup Options
opts = optiset('solver','mkltrnls','display','iter');
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'ydata',ydata,'ndec',3,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)

plot(Opt)

%% NLS2 w Weights
clc
%Function
i = (1:40)';
fun = @(x) x(1)*exp(-x(2)*i) + x(3);
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];
%Weighting
wts = ones(size(ydata)); wts([1 end]) = 1e3;
%Build & Solve
x0=[1.0; 0.0; 0.0];
Opt = opti('fun',fun,'ydata',ydata,'x0',x0);
OptW = opti('fun',fun,'ydata',ydata,'weights',wts,'x0',x0);
x = solve(Opt)
xw = solve(OptW)

plot(i,ydata,'ko',i(end),ydata(end),'ksq',i,fun(xw),'r*-',i,fun(x),'m*-')
xlim([length(i)*0.45 length(i)*1.01]); xlabel('i'); ylabel('y'); 
legend('Original Data','Target Point','Weighted Fit','Standard Fit');
title('NLS Curve Fit - Comparison of Weighted and Un-weighted');

%% NLS3 (Modified NLP - HS76)
clc
%Function
fun = @(x) [x(1);
            sqrt(0.5)*x(2);
            x(3);
            sqrt(0.5)*x(4);];
%Fitting Data
ydata = [0.0, 0.0, 0.0, 0.0];
%Constraints
lb = [0.0, 0.0, 0.0, 0.0];
ub = [inf, inf, inf, inf];
A = -[-1.0, -2.0, -1.0, -1.0;
     -3.0, -1.0, -2.0, 1.0];
b = -[-5.0, -0.4]';
Aeq = [0.0, 1.0, 4.0, 0.0];
beq = 1.5;
%Setup Options
opts = optiset('solver','levmar','display','iter');
%Build & Solve
x0 = [0.5, 0.5, 0.5, 0.5];
Opt = opti('fun',fun,'ydata',ydata,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'options',opts)
[x,fval,exitflag,info] = solve(Opt,x0)

%% NLS Example Prob
clc
%Get Problem
prob = nls_prob(19);
opts = optiset('display','iter');
%Build OPTI Object
Opt = opti(prob,opts);
%Solve
[x,fval,exitflag,info] = solve(Opt)
%Plot
plot(Opt,[],1)