function test_conf_formal()
% Fit Confidence & Statistic Testing
clc
clear

if (~exist('NonLinearModel.m','file'))
    fprintf('The statistics toolbox must be installed, skipping test\n');
    return;
end

%% Himmelblau Example
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','lmder','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
stats=fitStats(Opt,0.95,1);
stats.x = x;
% Repeat with Statistics Toolbox
mdl = NonLinearModel.fit(p,r,Rxn_rate,x0)   
% Check Solution
compareStats(mdl,stats,'Himmelblau');


%% Himmelblau Example with Weights
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
wts = ones(size(p)); wts(end) = 10;
%Setup Options
opts = optiset('solver','lmder','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'weights',wts,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
stats=fitStats(Opt,0.95,1);
stats.x = x;
% Repeat with Statistics Toolbox
mdl = NonLinearModel.fit(p,r,Rxn_rate,x0,'weights',wts)   
% Check Solution
compareStats(mdl,stats,'Himmelblau with Weights');


%% SAS Example 
% http://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_nlin_sect005.htm
clc
%Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
%Fitting Data
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [150 0.01];
Opt = opti('fun',eKin,'data',n,r,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
stats=fitStats(Opt,0.95,1);
stats.x = x;
% Repeat with Statistics Toolbox
mdl = NonLinearModel.fit(n,r,eKin,x0)
% Check Solution
compareStats(mdl,stats,'SAS Example');


% %% Stats Toolbox Example (OPTI gets a different solution...)
% clc
% load carbig
% X = [Horsepower,Weight];
% y = MPG;
% idx = ~isnan(X(:,1)) & ~isnan(y); X = X(idx,:); y = y(idx);
% modelfun = @(b,x)b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5);
% beta0 = [-50 500 -1 500 -1];
% %Build & Solve
% Opt = opti('fun',modelfun,'data',X,y,'x0',beta0,'opts',optiset('solver','nl2sol','display','off'));
% [x,fval,exitflag,info] = solve(Opt);
% stats=fitStats(Opt,0.95,1);
% stats.x = x;
% % Repeat with Statistics Toolbox
% mdl = NonLinearModel.fit(X,y,modelfun,beta0)
% % Check Solution
% compareStats(mdl,stats,'Stats Toolbox Example');



function compareStats(mdl,stats,str)

tol = 1e-5;
%Check Standard Stats
opticheckval.relErrorCheck(stats.Rsquare,mdl.Rsquared.Ordinary,sprintf('%s: R^2',str),tol);
opticheckval.relErrorCheck(stats.AdjRsquare,mdl.Rsquared.Adjusted,sprintf('%s: Adj R^2',str),tol);
opticheckval.relErrorCheck(stats.RMSE,mdl.RMSE,sprintf('%s: RMSE',str),tol);
opticheckval.relErrorCheck(stats.DFE,mdl.DFE,sprintf('%s: DFE',str),tol); 

%Parameter Stats
for i = 1:length(stats.Param.StdError)
    opticheckval.relErrorCheck(stats.Param.StdError(i),mdl.Coefficients.SE(i),sprintf('%s: SE[%d]',str,i),tol);
    opticheckval.relErrorCheck(stats.Param.tStat(i),mdl.Coefficients.tStat(i),sprintf('%s: tStat[%d]',str,i),tol);
end

%Confidence Interval
cis = coefCI(mdl); cis = cis(:,2) - mdl.Coefficients.Estimate; 
cio = stats.ConfInt;
for i = 1:length(stats.Param.StdError)
    opticheckval.relErrorCheck(cio(i),cis(i),sprintf('%s: CI[%d]',str,i),tol);
end

%Prediction Bounds
[~,ys] = predict(mdl,stats.ConfBnds.xdata,'Simultaneous','true');
yo = stats.ConfBnds.bnds;
idx = randi(numel(ys),min(10,numel(ys)));
for i = 1:length(idx)
    opticheckval.relErrorCheck(yo(idx(i)),ys(idx(i)),sprintf('%s: Pred[%d]',str,idx(i)),tol);
end