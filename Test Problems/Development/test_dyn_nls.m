%% Testing Dynamic Parameter Fitting
clear
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

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 2 Param, 1 State
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1) + p(2);
z0 = 1; %ic

%Generate measurement data
p = [2.5; 5.5];             %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement timess
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [1;1]; %inital parameter guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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

%Build OPTI Object
p0 = 1; %inital parameter guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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
zm = zm(:,state);

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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
zm = zm(:,state);

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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

tm = {tm1;tm2};
zm = {zm1;zm2};

%Build OPTI Object
theta0 = 1; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

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

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','matlab','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

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

tm = {tm1;tm2};
zm = {zm1;zm2};

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','matlab','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

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
zm = zm(:,state);

%Build OPTI Object
p0 = [1;0.1]; %inital parameter guess
dopts = optidynset('stateIndex',state,'sensitivity','nd');
opts = optiset('solver','lmder','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System
clc
%ODE to fit
ode = @(t,z,p) [p(1)*(z(2) - z(1));
                z(1)*(p(2) - z(3)) - z(2);
                z(1)*z(2) - p(3)*z(3)];
z0 = [5.7;10.50;30.58]; %ic

%Generate measurement data
p = [10,46,8/3];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:4;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('sensitivity','cs');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% Lorenz System [Analytical Derivatives]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*(z(2) - z(1));
                z(1)*(p(2) - z(3)) - z(2);
                z(1)*z(2) - p(3)*z(3)];
z0 = [5.7;10.50;30.58]; %ic
%Analytical Derivatives
dfdz = @(t,z,p) sparse([-p(1), p(1), 0;
                 p(2) - z(3), -1, -z(1);
                 z(2), z(1), -p(3)]);
dfdp = @(t,z,p) [z(2) - z(1), 0, 0;
                 0, z(1), 0;
                 0, 0, -z(3)];

%Generate measurement data
p = [10,46,8/3];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:4;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [Analytical Derivatives + Don't Measure 1st State]
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
tm = 0:0.1:4;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%States to Measure
state = [2 3];
zm = zm(:,state);

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',state,'sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [Analytical Derivatives + Non-Zero Starting Measurement Time]
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
tm = 0.5:0.1:4;             %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [Analytical Derivatives + Different Measurement Times + Non-Zero Initial Time (expected warning)]
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

tm = {tm1;tm2;tm3};
zm = {zm1;zm2;zm3};

%Given states 2 + 3 start from non-zero time, we should estimate
% z0(2:3) = NaN;,10.4,30.55

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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

tm = {tm1;tm2;tm3};
zm = {zm1;zm2;zm3};

%Given states 2 + 3 start from non-zero time, we should estimate
z0(2:3) = NaN;

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3,10.4,30.55]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [Analytical Derivatives + Different Measurement Times + State 2 and 3 + Non-Zero Initial Time + Estimate IC]
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
tm2 = 0.2:0.1:4;            %measurement times
tm3 = 0.5:0.15:3;              
[~,zm] = ode45(oi,[0 tm2],z0); zm2 = zm(2:end,2); %measurements
[~,zm] = ode45(oi,[0 tm3],z0); zm3 = zm(2:end,3); %measurements

tm = {tm2;tm3};
zm = {zm2;zm3};

%Estaimte Start 3 Initial Condition
z0(2:3) = NaN;

%States of Interest
stateIndex = [2 3];

%Build OPTI Object
p0 = [10+0.1;46+0.2;8/3;10.5;30.56]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'stateIndex',stateIndex);
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [INCORRECT Analytical Derivatives]
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
                 0, 1, -z(3)];

%Generate measurement data
p = [10,46,8/3];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:4;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('solver','auto','derivCheck','on','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Lorenz System [NO Sensitivity Integration]
clc
%ODE to fit
ode = @(t,z,p) [p(1)*(z(2) - z(1));
                z(1)*(p(2) - z(3)) - z(2);
                z(1)*z(2) - p(3)*z(3)];
z0 = [5.7;10.50;30.58]; %ic

%Generate measurement data
p = [10,46,8/3];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:4;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Build OPTI Object
p0 = [10+0.1,46+0.2,8/3]; %inital parameter guess
dopts = optidynset('sensitivity','none');
opts = optiset('solver','auto','derivCheck','on','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Flame ODE System
clc
%ODE to fit
ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;
z0 = 0.01; %ic

%Generate measurement data
p = 1.5;                    %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:1:2/z0;              %measurement times
[~,zm] = ode15s(oi,tm,z0);  %measurements

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('integrator','ode15s');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% Flame ODE System [Tricky DIW Version]
% clc
% %ODE to fit
% ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;
% z0 = 0.001; %ic
% 
% %Generate measurement data
% p = 1.5;                    %true parameter value
% oi = @(t,z) ode(t,z,p);     %ode integrator function
% tm = [0 65:70 80:0.1:90 100]*10;         %measurement times
% [~,zm] = ode15s(oi,tm,z0);  %measurements
% 
% %Build OPTI Object
% p0 = 1; %inital parameter guess
% dopts = optidynset('integrator','ode45');
% opts = optiset('display','iter','iterfun',@optiplotlogfval,'dynamicOpts',dopts);
% Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)
% 
% [x,f,e,i] = solve(Opt)
% plot(Opt)
% 
% %% Flame ODE System [Tricky DIW Version NO SENSITIVITY]
% clc
% %ODE to fit
% ode = @(t,z,p) p(1)*z(1)^2 - z(1)^3;
% z0 = 0.001; %ic
% 
% %Generate measurement data
% p = 1.5;                    %true parameter value
% oi = @(t,z) ode(t,z,p);     %ode integrator function
% tm = [0 65:70 80:0.1:90 100]*10;         %measurement times
% [~,zm] = ode15s(oi,tm,z0);  %measurements
% 
% %Build OPTI Object
% p0 = 1; %inital parameter guess
% dopts = optidynset('integrator','ode15s','sensitivity','nd');
% opts = optiset('display','iter','iterfun',@optiplotlogfval,'dynamicOpts',dopts);
% Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)
% 
% [x,f,e,i] = solve(Opt)
% plot(Opt)


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

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

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

%Estimate initial condition
z0 = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess
dopts = optidynset('sensitivity','nd');
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',theta0,'z0',z0,'options',opts)

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

%Replace z0 with NaN to indicate to estimate it
z0 = NaN;

%Build OPTI Object
theta0 = [1;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z01 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
theta0 = [1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z01 with NaN to indicate to estimate it
z0(2) = NaN;

%Build OPTI Object
theta0 = [1;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [1;0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% 0 Param, 2 State, + Solve for Both z0
clc
%ODE to fit
ode = @(t,z,p) [z(1) + z(2); z(1)];
z0 = [1;3]; %ic

%Generate measurement data
p = 0;                  %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.2:1;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',theta0,'z0',z0,'options',opts)

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

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
theta0 = [1;0.5;2.5]; %inital parameter guess + initial state guess
dopts = optidynset('integrator','ode45','sensitivity','nd');
opts = optiset('display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z0 with NaN to indicate to estimate it
z0(1) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z0 with NaN to indicate to estimate it
z0(2) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.1;0.1;0.1]; %inital parameter guess + initial state guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
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
opts = optiset('solver','ipopt','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

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
plot(Opt)


%% DIW MAPLE ODEs [IC Estimate]
clc
ode = @(t,z,p) [-p(1)*z(1) + 4; 
                2*z(1) - p(1)*z(2) + 5; 
                -4*z(1) - 2*z(3) - p(2)];
z0 = [-1.5;1.25;1]; %ic

%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:2;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0([1, 3]) = NaN;

%Build OPTI Object
theta0 = [1;0.1;1;0.1]; %inital parameter guess + initial state guess
dopts = optidynset('sensitivity','none');
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:2;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + States 1 3]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:2;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%States of Interest
states = [1 3];
zm = zm(:,states);

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',states);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + States 1 3 + Repeated Measurements]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:2;               %measurement times
[~,zm1] = ode45(oi,tm,z0);   %measurements
[~,zm2] = ode45(oi,tm,z0*0.9);   %measurements

%Concatenate measurements
tm_m = [tm tm]';
zm_m = [zm1;zm2];

%States of Interest
states = [1 3];
zm_m = zm_m(:,states);

%Replace z0 with NaN to indicate to estimate it
z0(1:3) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',states);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_m,zm_m,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + States 1 3 + Different Measurement Times]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0:0.1:2;              %measurement times
tm3 = 0:0.05:2;             
[~,zm] = ode45(oi,tm1,z0); zm1 = zm(:,1);  %measurements
[~,zm] = ode45(oi,tm3,z0); zm3 = zm(:,3);  

tm = {tm1,tm3};
zm = {zm1,zm3};

%States of Interest
states = [1 3];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',states);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + States 1 3 + Different Measurement Times + Repeated Points]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0:1:2;              %measurement times
tm3 = 0:0.5:2;             
[~,zm] = ode45(oi,tm1,z0); zm1 = zm(:,1);  %measurements
[~,zm] = ode45(oi,tm3,z0); zm3 = zm(:,3);  
[~,zm] = ode45(oi,tm3,z0*0.9); zm3b = zm(:,3);  

tm = {tm1,[tm3 tm3 ]};
zm = {zm1,[zm3;zm3b]};

%States of Interest
states = [1 3];

%Replace z0 with NaN to indicate to estimate it
z0([1 3]) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',states,'sensitivity','user');
opts = optiset('solver','lmder','display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + States 1 3 + Different Measurement Times + Non-Zero Initial Time]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm1 = 0.5:0.5:2;              %measurement times
tm3 = 0:0.05:2;             
[~,zm] = ode45(oi,[0 tm1],z0); zm1 = zm(2:end,1);  %measurements
[~,zm] = ode45(oi,tm3,z0); zm3 = zm(:,3);  

tm = {tm1,tm3};
zm = {zm1,zm3};

%States of Interest
states = [1 3];

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'stateIndex',states);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'options',opts)

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
[~,zm1b] = ode45(odeInt,[0 tm1],z0*0.9); zm1b = zm1b((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
 tm_multi = {[tm1 tm1];tm2;tm3};
 zm_multi = {[zm1;zm1b];zm2;zm3};

% We will need to estimate all initial conditions in this problem
 z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

%Build OPTI Object
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'sensitivity','user');
opts = optiset('solver','matlab','display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + Different Measurement Times + Repeated Points [2x states] + Non-Zero Initial Time]
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
             
p = [2.345;1.1];            %true parameter value
% Measurement Times for each State
 tm1  = 0.5:0.1:2;              %state 1
 tm2  = 1:0.1:2;                %state 2
 tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
odeInt = @(t,z) ode(t,z,p);     %ode integrator function
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm1b] = ode45(odeInt,[0 tm1],z0*0.4); zm1b = zm1b((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm2b] = ode45(odeInt,[0 tm2],z0*0.3); zm2b = zm2b((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);

% Group measurements and time stamps in cell arrays
 tm_multi = {[tm1 tm1];[tm2 tm2];tm3};
 zm_multi = {[zm1;zm1b];[zm2;zm2b];zm3};

% We will need to estimate all initial conditions in this problem
 z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

%Build OPTI Object
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'sensitivity','nd');
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + Different Measurement Times + Repeated Points [3x states] + Non-Zero Initial Time]
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

p = [2.345;1.1];            %true parameter value
% Measurement Times for each State
 tm1  = 0.5:0.1:2;              %state 1
 tm2  = 1:0.1:2;                %state 2
 tm3  = 0.2:0.2:2;              %state 3

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
odeInt = @(t,z) ode(t,z,p);     %ode integrator function
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm1b] = ode45(odeInt,[0 tm1],z0*0.4); zm1b = zm1b((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm2b] = ode45(odeInt,[0 tm2],z0*0.3); zm2b = zm2b((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);
[~,zm3b] = ode45(odeInt,[0 tm3],z0*0.95); zm3b = zm3b((2:end),3);
[~,zm3c] = ode45(odeInt,[0 tm3],z0*1.1); zm3c = zm3c((2:end),3);

% Group measurements and time stamps in cell arrays
 tm_multi = {[tm1 tm1];[tm2 tm2];[tm3 tm3 tm3]};
 zm_multi = {[zm1;zm1b];[zm2;zm2b];[zm3;zm3b;zm3c]};

% We will need to estimate all initial conditions in this problem
 z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

%Build OPTI Object
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'sensitivity','nd');
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW MAPLE ODEs [Analytical Derivatives + IC Estimate + Different Measurement Times + Repeated Points [3x states] + Non-Zero Initial Time]
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

p = [2.345;1.1];            %true parameter value             
% Measurement Times for each State
 tm1  = 0.5:0.1:2;              %state 1
 tm2  = 1:0.1:2;                %state 2
 tm2b  = 1:0.03:2;             
 tm3  = 0.2:0.2:2;              %state 3
 tm3c = 0.8:0.2:2;  

% Solve ODEs and index Measurement Data
% Note 0 is required below just to match initial condition in this example
odeInt = @(t,z) ode(t,z,p);     %ode integrator function
[~,zm1] = ode45(odeInt,[0 tm1],z0); zm1 = zm1((2:end),1);
[~,zm1b] = ode45(odeInt,[0 tm1],z0*0.4); zm1b = zm1b((2:end),1);
[~,zm2] = ode45(odeInt,[0 tm2],z0); zm2 = zm2((2:end),2);
[~,zm2b] = ode45(odeInt,[0 tm2b],z0*0.8); zm2b = zm2b((2:end),2);
[~,zm3] = ode45(odeInt,[0 tm3],z0); zm3 = zm3((2:end),3);
[~,zm3b] = ode45(odeInt,[0 tm3],z0*0.95); zm3b = zm3b((2:end),3);
[~,zm3c] = ode45(odeInt,[0 tm3c],z0*1.1); zm3c = zm3c((2:end),3);

% Group measurements and time stamps in cell arrays
 tm_multi = {[tm1 tm1];[tm2 tm2b];[tm3 tm3 tm3c]};
 zm_multi = {[zm1;zm1b];[zm2;zm2b];[zm3;zm3b;zm3c]};

% We will need to estimate all initial conditions in this problem
 z0 = [NaN;NaN;NaN];

% New Initial Guess Vector [p;z0];
 theta0 = [1;0.5;0.5;0.5;0.5];

%Build OPTI Object
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'initialT',0,'sensitivity','user');
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm_multi,zm_multi,'x0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% DIW MAPLE ODEs [Analytical Derivatives + Stiff Solver]
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
             
%Generate measurement data
p = [2.345;1.1];            %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tm = 0:0.1:2;               %measurement times
[~,zm] = ode45(oi,tm,z0);   %measurements

%Replace z0 with NaN to indicate to estimate it
z0(1:2) = NaN;

%Build OPTI Object
p0 = [1;0.5;1;0.5]; %inital parameter guess + initial state guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp);
opts = optiset('display','iter','derivCheck','on','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm(:),'x0',p0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% DIW Problem File
clc
% Duffing's equation 
a = 0.5; b = 1; c = 0.05; % CONSTANT parameters 
A = 3; % Adjustable parameter

%ODE to fit
ode = @(t,z,p) [z(2); ... 
                -c*z(2) - p(3)*z(1) - p(2)*z(1)^3 + p(1)*cos(t)]; 
z0 = [1; -2]; %ic

%Analytical Derivatives
dfdz = @(t,z,p) [0,                 1;
                 -p(3)-3*p(2)*z(1)^2, -c];
dfdp = @(t,z,p) [0,      0,        0; 
                 cos(t), -z(1)^3,  -z(1)];

%Generate measurement data
p = [A; b; a];             %true parameter value
oi = @(t,z) ode(t,z,p);     %ode integrator function
tspan = [0,40];               %measurement times
[tm,zm] = ode45(oi,tspan,z0);   %measurements

n = length(tm); 
N = min(100,n); 
idx = randi(n,N,2); 
% i1 = sort(idx(:,1)); i2 = sort(idx(:,2)); 
i1 = unique(idx(:,1)); i2 = unique(idx(:,2));
 
 
tm1 = tm(i1,1); zm1 = zm(i1,1); 
tm2 = tm(i2); zm2 = zm(i2,2); 
 
tm = {tm1,tm2}; zm = {zm1,zm2};

%plot(tm, zm,'o-')  % raw data 
 
z0est = NaN*z0;
theta = [p;z0]
theta0 = theta*0.95; %inital parameter guess
dopts = optidynset('dfdz',dfdz,'dfdp',dfdp,'integrator','ode15s','sensitivity','user');
opts = optiset('solver','auto','display','iter','derivCheck','on','dynamicOpts',dopts);

Opt = opti('ode',ode,'data',tm,zm,'theta0',theta0,'z0',z0est,'options',opts)

tic; 
[param_fit,f,e,i] = solve(Opt)
toc
%plot(Opt)

%
        
ofit = @(t,z) ode(t,z,param_fit);     %ode integrator function
          
[tfit,zfit] = ode45(ofit,[tspan(1), tspan(end)],z0);   %measurements
ofit = @(t,z) ode(t,z,theta0);  
[tfit0,zfit0] = ode45(ofit,[tspan(1), tspan(end)],theta0(end-1:end));   %measurements

%
subplot(2,1,1); 
h = plot(tm1, zm1,'o', tfit, zfit(:,1),'-', tfit0, zfit0(:,1),'-');  % raw data 
set(h(3),'color',0.8*[1,1,1]); 
errf = norm(theta-param_fit); 
title(sprintf('A true = %2g, Afit = %2g, b true = %2.2g, bfit = %2.2g, Norm(err) = %2.2f %%',A,param_fit(1), b, param_fit(2), errf))
% set(gca,'XTickl',''); 
 
subplot(2,1,2); 
h = plot(tm2, zm2,'o', tfit, zfit(:,2),'-', tfit0, zfit0(:,2),'-');  % raw data 
set(h(3),'color',0.8*[1,1,1]); 


%% Analyical Derivatives1
% ode = @(t,z,p) [p(1)*(z(2) - z(1));
%                 z(1)*(p(2) - z(3)) - z(2);
%                 z(1)*z(2) - p(3)*z(3)];
syms p1 p2 p3 z1 z2 z3
o = [p1*(z2 - z1); z1*(p2 - z3) - z2; z1*z2 - p3*z3]

jacobian(o,[z1 z2 z3])
jacobian(o,[p1 p2 p3])

%% Analyical Derivatives2
ode = @(t,z,p) [-p(1)*z(1) + 4; 
                2*z(1) - p(1)*z(2) + 5; 
                -4*z(1) - 2*z(3) - p(2)];
syms p1 p2 z1 z2 z3
o = [-p1*z1 + 4; 2*z1 - p1*z2 + 5; -4*z1 - 2*z3 - p2];

jacobian(o,[z1 z2 z3])
jacobian(o,[p1 p2])

%% Symbolic Testing
clc
ode = @(t,z,p) [z(2); 
                p(1)*(1-z(1)^2)*z(2) - z(1)];
   
z0 = [1.5,2]'; % Initial Conditions
p = 3;     % True Parameter Values, mu            

%Generate Partial Derivatives
[dfdz,dfdp] = symDynJac(ode)            

%Comparison
dfdz(1,z0,p)
dfdp(1,z0,p)

mklJac(@(z) ode(1,z,p),z0)
mklJac(@(p) ode(1,z0,p),p)

%% Symbolic Testing 2
clc
a = 0.5; b = 1; c = 0.05; % CONSTANT parameters 
A = 3; % Adjustable parameter
%ODE to fit
ode = @(t,z,p) [z(2); ... 
                -c*z(2) - p(3)*z(1) - p(2)*z(1)^3 + p(1)*cos(t)]; 
z0 = [1; -2]; %ic
p = [A; b; a];             %true parameter value

%Generate Partial Derivatives
[dfdz,dfdp] = symDynJac(ode)            

%Comparison
dfdz(1,z0,p)
dfdp(1,z0,p)

mklJac(@(z) ode(1,z,p),z0)
mklJac(@(p) ode(1,z0,p),p)





