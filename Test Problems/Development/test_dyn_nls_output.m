%% Testing Dynamic Parameter Fitting with output function
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

% Analytic derivs
dfdz = @(t,z,p) p(1);
dfdp = @(t,z,p) z(1);

%Build OPTI Object
p0 = 1; %inital parameter guess
dopts = optidynset('dfdp',dfdp,'dfdz',dfdz);
opts = optiset('solver','auto','display','iter','dynamicOpts',dopts);
Opt = opti('ode',ode,'data',tm,zm,'x0',p0,'z0',z0,'bounds',0,2,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% Gradient comparison
clc
f = Opt.nlprob.fun
g = Opt.nlprob.grad

% Output
[~,zm] = ode45(oi,tm,z0);
h = @(z,p) 2*z + 3*p;

h(


function 