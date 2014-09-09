%% 1 Param, 1 State
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(end) = zm(end)*1.1;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 1 Param, 1 State, End Point Weighted [No Sensitivity]
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(end) = zm(end)*1.1;

weights = ones(size(zm)); weights(end) = 100;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','none');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 1 Param, 1 State, End Point Weighted [With Sensitivity]
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(end) = zm(end)*1.1;

weights = ones(size(zm)); weights(end) = 100;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 1 Param, 1 State, Reversed Timestamps
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm = flipud(zm); tm = fliplr(tm);

zm(3) = zm(3)*1.5;
weights = ones(size(zm)); weights(3) = 0;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% 1 Param, 1 State, End Point Weighted [No Sensitivity] as NLP
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(end) = zm(end)*1.1;

weights = ones(size(zm)); weights(end) = 100;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','none');
opts = optiset('solver','ipopt','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 1 Param, 1 State, End Point Weighted [With Sensitivity] as NLP
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1);
z0 = 1; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements
zm(end) = zm(end)*1.1;

weights = ones(size(zm)); weights(end) = 100;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','ipopt','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


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
zm(4,1) = zm(4,1)*1.5;

weights = ones(size(zm)); weights(4,1) = 0;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

%State to Measure
state = 1;
zm = zm(:,state); zm(4) = zm(4)*1.5;

weights = ones(size(zm)); weights(4) = 0;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

%State to Measure
state = 2;
zm = zm(:,state); zm(4) = zm(4)*1.5;

weights = ones(size(zm)); weights(4) = 0;

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

zm1(4) = zm1(4)*1.5;

tm = {tm1;tm2};
zm = {zm1;zm2};

wts1 = ones(size(zm1)); wts1(4) = 0;
wts2 = ones(size(zm2));
weights = {wts1;wts2};

%Build OPTI Object
theta0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

tm = fliplr(tm);
zm = flipud(zm);
zm(3,1) = zm(3,1)*1.5;

weights = ones(size(zm));
weights(3,1) = 0;

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

zm2(4) = zm2(4)*0.1;
tm = {tm1;tm2};
zm = {zm1;zm2};

wts1 = ones(6,1);
wts2 = ones(size(zm2)); wts2(4) = 0;
weights = {wts1;wts2};

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

%State to Measure
state = 2;
zm = zm(:,state); zm(4) = zm(4)*1.5;

weights = ones(size(zm)); weights(4) = 0;

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','lmder','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% Lorenz System [Analytical Derivatives + Different Measurement Times + Non-Zero Initial Time + Estimate IC]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*(z(2) - z(1));
                z(1)*(p(2) - z(3)) - z(2);
                z(1)*z(2) - p(3)*z(3)];
z0 = [5.7;10.50;30.58]; %ic
%Analytical Derivatives
dfdz = @(t,z,p) [-p(1), p(1), 0;
                 p(2) - z(3), -1, -z(1);
                 z(2), z(1), -p(3)];
dfdp = @(t,z,p) [z(2) - z(1), 0, 0;
                 0, z(1), 0;
                 0, 0, -z(3)];

%Generate measurement data
p = [10,46,8/3];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0:0.2:2;
tm2 = 0.2:0.1:4;            %measurement times
tm3 = 0.5:0.5:3;              
[~,zm] = ode45(oi,tm1,z0); zm1 = zm(:,1); %measurements
[~,zm] = ode45(oi,[0 tm2],z0); zm2 = zm(2:end,2); %measurements
[~,zm] = ode45(oi,[0 tm3],z0); zm3 = zm(2:end,3); %measurements

zm3(3) = zm3(3)*1.2;

tm = {tm1;tm2;tm3};
zm = {zm1;zm2;zm3};

wts1 = ones(size(zm1));
wts2 = ones(size(zm2));
wts3 = ones(size(zm3)); wts3(3) = 0;
weights = {wts1;wts2;wts3};

%Given states 2 + 3 start from non-zero time, we should estimate
z0(2:3) = NaN;

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3,10.4,30.55]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('solver','nl2sol','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

wts1 = 2*ones(size(zm1));
wts2 = ones(size(zm2)); 
weights = [wts1;wts2];

%Estimate initial condition
z0 = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','nl2sol','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

wts1 = 2*ones(size(zm1));
wts2 = ones(size(zm2)); 
weights = [wts1;wts2];

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

wts1 = 2*ones(size(zm1));
wts2 = ones(size(zm2)); 
weights = [wts1;wts2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [1;0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

wts1 = 2*ones(size(zm1));
wts2 = ones(size(zm2)); 
weights = [wts1;wts2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%States to Measure
state = 1;
zm_m = zm_m(:,state);
weights = weights(:,state);

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
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

wts1 = 2*ones(size(zm1));
wts2 = ones(size(zm2)); 
weights = [wts1;wts2];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%States to Measure
state = 2;
zm_m = zm_m(:,state);
weights = weights(:,state);

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + Different Measurement Times + Repeated Points + Non-Zero Initial Time]
clc
ode = @(t,z,p) [-p(1)*z(1) + 4; 
                2*z(1) - p(1)*z(2) + 5; 
                -4*z(1) - 2*z(3) - p(2)];
z0 = [-1.5;1.25;1]; %ic

%Analytical Derivatives
dfdz = @(t,z,p) [-p(1), 0, 0;
                 2, -p(1), 0;
                 -4, 0, -2];
dfdp = @(t,z,p) [-z(1), 0
                 -z(2), 0
                 0, -1];
             
% Measurement Times for each State
 tm1  = 0.5:0.1:2;              %state 1
 tm2  = 1:0.1:2;                %state 2
 tm3  = 0.2:0.2:2;              %state 3

 p = [2.345;1.1];            %true parameter value
% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
odeInt = @(t,z) ode(t,z,p);     %ode integrator function
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm1b] = ode45(odeInt,[0 tm1],z0*0.5); zm1b = zm1b((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
 tm_multi = {[tm1 tm1];tm2;tm3};
 zm_multi = {[zm1;zm1b];zm2;zm3};
 
 wts1 = ones(size(zm1));
 wts1b = 10*ones(size(zm1));
 wts2 = ones(size(zm2)); 
 wts3 = ones(size(zm3)); 
 weights = {[wts1;wts1b];wts2;wts3};

% We will need to estimate all initial conditions in this problem
 z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

%Build OPTI Object
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'sensitivity','user');
opts = optiset('solver','auto','display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'x0',theta0,'z0',z0,'weights',weights,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


