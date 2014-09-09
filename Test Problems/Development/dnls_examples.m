%% DNLS Example1
clc
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


%% DNLS Example2
clc
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

%% Example 3
clc

% Add Bounds
 lb = [1.5;0;0.1;0.1];
 ub = [3.5;2;1.5;1.5];

% Create OPTI Object
 Opt = opti('ode',ode,'data',tm,zm,'bounds',lb,ub,'z0',z0,'theta0',theta0)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)

%% DNLS Example4
clc
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

%% DNLS Example5
clc
% Initial Conditions
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

% Indicate Initial Conditions to Solve for Using NaNs
 z0 = [-1.5;NaN;NaN];

% Create OPTI Object
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)

%% DNLS Example6
clc
% Initial Conditions
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm1  = 0.5:0.1:2;   %state 1
tm2  = 1:0.1:2;    %state 2
tm3  = 0.2:0.2:2;   %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
tm_multi = {tm1;tm2;tm3};
zm_multi = {zm1;zm2;zm3};

% Indicate Initial Conditions to Solve for Using NaNs
 z0 = [NaN;NaN;NaN];
 
% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

% Create OPTI Object
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)


%% DNLS Example7
clc
% Initial Conditions
z0 = [-1.5;1.25;1];

% Measurement Times for each State
tm1  = 0.5:0.1:2;   %state 1
tm1b = 1:0.25:2;   %state 1
tm2  = 1:0.1:2;    %state 2
tm3  = 0.2:0.2:2;   %state 3

%ODE with modified parameters (representing a different run)
odeIntB = @(t,z) ode(t,z,p*0.85);  

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
[~,zm1]  = ode45(odeInt,[0 tm1],z0);  zm1  = zm1((2:end),1);
[~,zm1b] = ode45(odeIntB,[0 tm1b],z0); zm1b = zm1b((2:end),1);
[~,zm2]  = ode45(odeInt,[0 tm2],z0);  zm2  = zm2((2:end),2);
[~,zm3]  = ode45(odeInt,[0 tm3],z0);  zm3  = zm3((2:end),3);
[~,zm3b] = ode45(odeInt,[0 tm3],z0*0.75); zm3b = zm3b((2:end),3);

% Group measurements and time stamps in cell arrays
tm_multi = {[tm1 tm1b];tm2;[tm3 tm3]};
zm_multi = {[zm1;zm1b];zm2;[zm3;zm3b]};

% Indicate Initial Conditions to Solve for Using NaNs
 z0 = [NaN;NaN;NaN];
 
% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

% Create OPTI Object
dopts = optidynset('initialT',0);
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)


%% DNLS Example8
clc
% Declare Symbolic Variables
syms p1 p2 z1 z2 z3
% Symbolic ODE
ODE = [-p1*z1 + 4; 
      2*z1 - p1*z2 + 5; 
      -4*z1 - 2*z2*z3 - p2];

% Evaluate Jacobians
dfdz_sym = jacobian(ODE,[z1 z2 z3])
dfdp_sym = jacobian(ODE,[p1 p2])

[dfdz,dfdp] = symDynJac(ode)


% Analytical Derivative Expressions
 dfdz = @(t,z,p) [-p(1) 0,       0;
                  2,    -p(1),   0;
                  -4,   -2*z(3), -2*z(2)];
 dfdp = @(t,z,p) [-z(1), 0;
                  -z(2), 0;
                  0,    -1];

% Set Dynamic Options
 dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'sensitivity','User','initialT',0);

% General OPTI Options
 opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);

% Create OPTI Object
  Opt = opti('ode',ode,'data',tm_multi,zm_multi,'z0',z0,'theta0',theta0,...
             'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)

% Plot the Solution
plot(Opt)


%% DNLS Ex9
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;
z0 = 0.001; %ic

%Generate measurement data
p = 1.5;                          %true parameter value
oi = @(t,z) ode(t,z,p);           %ode integrator function
tm = [0 65:70 80:1:90 100]*10;  %measurement times
[~,zm] = ode15s(oi,tm,z0);        %measurements

% Starting Guess
theta0 = 1;

%Build OPTI Object
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
subplot(121)
plot(OptNoDer); title('No Sensitivity');
subplot(122)
plot(OptWDer); title('With Sensitivity');

%% DNS Ex9b
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

% Plot
plot(Opt)


%% DNS Ex9c
clc
% Set Dynamic Options
dopts = optidynset('integrator','ode15s','sensitivity','None');

% General OPTI Options
opts = optiset('solver','nomad','display','iter','dynamicOpts',dopts);

% Create OPTI Object
Opt = opti('ode',ode,'data',tm,zm,'bounds',0.5,2.5,'z0',z0,'theta0',theta0,'options',opts)

% Solve
[theta,fval,exitflag,info] = solve(Opt)