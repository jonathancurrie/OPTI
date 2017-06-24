%% Confidence Limit Testing
% J.Currie May 2014
clc
clear
% format compact
% format short g

haveStatsToolbox = exist('NonlinearModel.m','file');

%% Himmelblau Example
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
s=fitStats(Opt,0.95,1);
%Plot
clf
plot(Opt)

% Repeat with Statistics Toolbox
if (haveStatsToolbox)
    mdl = NonLinearModel.fit(p,r,Rxn_rate,x0)   
end

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData] = prepareCurveData( p, r );
    % Set up fittype and options.
    ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end

%% Himmelblau Example with weights
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
wts = ones(size(p)); wts(end) = 10;
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'weights',wts,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
stats = fitStats(Opt,0.95,1);
%Plot
plot(Opt)

% Repeat with Statistics Toolbox
if (haveStatsToolbox)
    mdl = NonLinearModel.fit(p,r,Rxn_rate,x0,'weights',wts)
end

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData, weights] = prepareCurveData( p, r, wts );
    % Set up fittype and options.
    ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    opts.Weights = weights;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end


%% Himmelblau Example WITH BOUNDS [LB]
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'bounds',[0;0],[5e-3,1],'options',opts);
[x,fval,exitflag,info] = solve(Opt);
fitStats(Opt,0.95);
%Plot
plot(Opt)

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData] = prepareCurveData( p, r );
    % Set up fittype and options.
    ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end

%% Himmelblau Example WITH BOUNDS [UB]
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'bounds',[0;0.03],[7e-3,1],'options',opts);
[x,fval,exitflag,info] = solve(Opt);
fitStats(Opt,0.95);
%Plot
plot(Opt)

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData] = prepareCurveData( p, r );
    % Set up fittype and options.
    ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end

%% Himmelblau Example WITH BOUNDS [LB,UB]
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'bounds',[0;0.03],[5e-3,1],'options',opts);
[x,fval,exitflag,info] = solve(Opt);
fitStats(Opt,0.95);
%Plot
plot(Opt)

%% Himmelblau Example WITH BOUNDS [UB] and WEIGHTS
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
wts = ones(size(p)); wts(end) = 100;
%Setup Options
opts = optiset('solver','auto','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'weights',wts,'x0',x0,'bounds',[0;0],[4e-3,1],'options',opts);
[x,fval,exitflag,info] = solve(Opt);
s=fitStats(Opt,0.95,1);
[x-s.ConfInt x+s.ConfInt]
%Plot
plot(Opt)

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData, weights] = prepareCurveData( p, r, wts );
    % Set up fittype and options.
    ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    opts.Weights = weights;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end


%% SAS Example 
% http://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_nlin_sect005.htm
clc
%Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
%Fitting Data
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','mkltrnls','display','off');
%Build & Solve
x0 = [150 0.01];
Opt = opti('fun',eKin,'data',n,r,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
s=fitStats(Opt,0.95,1);
%Plot
plot(Opt)

% Repeat with Statistics Toolbox
if (haveStatsToolbox)
    mdl = NonLinearModel.fit(n,r,eKin,x0)
end

%Compare with curve fitting toolbox if available
if(~isempty(which('cfit')))
    [xData, yData] = prepareCurveData( n, r );
    % Set up fittype and options.
    ft = fittype( 'a*x./(b+x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = Opt.prob.lb;
    opts.StartPoint = Opt.prob.x0;
    opts.TolFun = 1e-09;
    opts.TolX = 1e-09;
    opts.Upper = Opt.prob.ub;
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fitresult
    gof
    xd = s.ConfBnds.xdata;
    y = feval(fitresult,xd);
    ybnds = predint(fitresult,xd,0.95,'Functional','on');
    hold on
    plot(xd,y,':',xd,ybnds,':');
    hold off
end

%% Stats Toolbox Example
clc
if (haveStatsToolbox)
    load carbig
    X = [Horsepower,Weight];
    y = MPG;
    idx = ~isnan(X(:,1)) & ~isnan(y); X = X(idx,:); y = y(idx);
    modelfun = @(b,x)b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5);
    beta0 = [-50 500 -1 500 -1];

    %Build & Solve
    Opt = opti('fun',modelfun,'data',X,y,'x0',beta0,'opts',optiset('solver','auto','display','off'));
    [x,fval,exitflag,info] = solve(Opt);
    fitStats(Opt,0.95);

    % Repeat with Statistics Toolbox
    mdl = NonLinearModel.fit(X,y,modelfun,beta0)
end

%% Himmelblau Example Nonlinear Programming Solver
clc
%Function
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
%Fitting Data
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
%Setup Options
opts = optiset('solver','ipopt','display','off');
%Build & Solve
x0 = [5e-3 2e-2];
Opt = opti('fun',Rxn_rate,'data',p,r,'x0',x0,'options',opts);
[x,fval,exitflag,info] = solve(Opt);
fitStats(Opt,0.95);
%Plot
plot(Opt)

% Repeat with Statistics Toolbox
if (haveStatsToolbox)
    mdl = NonLinearModel.fit(p,r,Rxn_rate,x0)   
end



%% ODE1
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+.3*randn(size(zm));

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'bounds',0,2,'options',opts);

[x,f,e,i] = solve(Opt);
fitStats(Opt);
plot(Opt)

% Repeat with Statistics Toolbox
% mdl = fitnlm(tm,zm,,p0)

%% ODE1 [DIW Mod]
clc
%ODE to fit
ode = @(t,z,p) -p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = [0:0.1:3 5];               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(1:end-1) = zm(1:end-1)+.1*randn(size(zm(1:end-1)));

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'bounds',0,2,'options',opts);

[x,f,e,i] = solve(Opt);
fitStats(Opt);
plot(Opt)

%% NLS Example Prob
clc
%Get Problem
prob = nls_prob(19);
opts = optiset('display','iter','solver','lmder');
%Build OPTI Object
Opt = opti(prob,opts);
%Solve
[x,fval,exitflag,info] = solve(Opt);
fitStats(Opt,0.95);
%Plot
plot(Opt,[],1)

%% 1 Param, 2 State
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1.5*randn(size(zm));

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 2 State [only fitting first state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)]; %note z1 in 2nd ode allows us to estimate p1
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+.8*randn(size(zm));

%State to Measure
state = 1;
zm = zm(:,state);

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 2 State [only fitting second state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)]; %note z1 in 2nd ode allows us to estimate p1
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+.8*randn(size(zm));

%State to Measure
state = 2;
zm = zm(:,state);

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 2 State [different measurement times]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0:0.2:1;              %measurement times
tm2 = 0:0.15:1;
[~,zm] = ode45(oi,tm1,z0); zm1 = zm(:,1); %measurements
[~,zm] = ode45(oi,tm2,z0); zm2 = zm(:,2); %measurements

tm = {tm1;tm2};
zm = {zm1;zm2};

%Build OPTI Object
theta0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+.8*randn(size(zm));

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','matlab','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State [Different Measurement Times]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0:0.2:1;              %measurement times
tm2 = 0:0.15:1;
[~,zm] = ode45(oi,tm1,z0); zm1 = zm(:,1); %measurements
[~,zm] = ode45(oi,tm2,z0); zm2 = zm(:,2); %measurements

tm = {tm1;tm2};
zm = {zm1;zm2};

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','matlab','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State [Only fit 2nd state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+.8*randn(size(zm));

%State to Measure
state = 2;
zm = zm(:,state);

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','lmder','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 1 State, + Solve for z0
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+0.2*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 1 Param, 1 State [Repeated Measurements + z0]
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.95);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Estimate initial condition
z0 = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 1 State, + Solve for z0 + Perturb Initial Measurement (incorrect z0)
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;
%Perturb Initial zm
zm(1) = zm(1) + 0.5;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 1 State, + Solve for z0
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1) + p(2);
z0 = 1; %ic

%Generate measurement data
p = [2.5; 5.5];             %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+2*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 2 Param, 1 State [Repeated Measurements + z0]
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1) + p(2);
z0 = 1; %ic

%Generate measurement data
p = [2.5; 5.5];             %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.5);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 1 Param, 2 State, + Solve for z01
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z01 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 1 Param, 2 State, + Solve for z02
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z01 with NaN to indicate to estimate it
z0(2) = NaN;

%Build OPTI Object
theta0 = [1;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 2 State, + Solve for Both z0
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [1;0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 1 Param, 2 State, Repeated Measurements + Both z0
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [1;0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 1 Param, 2 State, Repeated Measurements + z01
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 1 Param, 2 State, Repeated Measurements + z02
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 2.345;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(2) = NaN;

%Build OPTI Object
theta0 = [1;1.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for z01
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for z01 + Repeated Measurements
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)

%% 2 Param, 2 State + Solve for z02
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0(2) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for both z0
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = zm+1*randn(size(zm));

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for both z0 + Repeated Measurements
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for both z0 [Only fit 2nd state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%States to Measure
state = 2;
zm = zm(:,state);

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for both z0 + Repeated Measurements [Only fit 1st state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%States to Measure
state = 1;
zm_m = zm_m(:,state);

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)


%% 2 Param, 2 State + Solve for both z0 + Repeated Measurements [Only fit 2nd state]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*z(1) + z(2); p(2)*z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm1] = ode45(oi,tm,z0);        %measurements RUN 1
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements RUN 2

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%States to Measure
state = 2;
zm_m = zm_m(:,state);

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
fitStats(Opt);
plot(Opt)
