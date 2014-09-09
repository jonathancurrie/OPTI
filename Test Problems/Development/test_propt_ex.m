%% TOMLAB PROPT Parameter Estimation Problems
% http://tomdyn.com/parameter_estimation_dynamic_systems.html
clc
clear

%% http://tomdyn.com/examples/catalyticCracking.html
clc
ode = @(t,z,p) [-(p(1)+p(3))*z(1).^2
                p(1)*z(1).^2-p(2)*z(2)];
z0 = [1;0]; %ic


% Various constants and expressions
y1meas = [1.0;0.8105;0.6208;0.5258;0.4345;0.3903;...
    0.3342;0.3034;0.2735;0.2405;0.2283;0.2071;0.1669;...
    0.153;0.1339;0.1265;0.12;0.099;0.087;0.077;0.069];
y2meas = [0;0.2;0.2886;0.301;0.3215;0.3123;0.2716;...
    0.2551;0.2258;0.1959;0.1789;0.1457;0.1198;0.0909...
    ;0.0719;0.0561;0.046;0.028;0.019;0.014;0.010];
tmeas = [0;0.025;0.05;0.075;0.1;0.125;...
    0.15;0.175;0.2;0.225;0.25;0.3;0.35;0.4;...
    0.45;0.5;0.55;0.65;0.75;0.85;0.95];

%Build OPTI Object
theta0 = [10;8;1]; %inital parameter guess
opts = optiset('display','iter');
Opt = opti('ode',ode,'data',tmeas,[y1meas' y2meas'],'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% http://tomdyn.com/examples/isometrizationAlpha.html
clc
ode = @(t,z,p) [-(p(1)+p(2))*z(1)
                p(1)*z(1)
                p(2)*z(1)-(p(3)+p(4))*z(3)+p(5)*z(5)
                p(3)*z(3)
                p(4)*z(3)-p(5)*z(5)];
z0 = [100;0;0;0;0]; %ic


y1meas = [88.35; 76.4; 65.1; 50.4; 37.5; 25.9; 14.0; 4.5];
y2meas = [7.3; 15.6; 23.1; 32.9; 42.7; 49.1; 57.4; 63.1];
y3meas = [2.3; 4.5; 5.3; 6.0; 6.0; 5.9; 5.1; 3.8];
y4meas = [0.4; 0.7; 1.1; 1.5; 1.9; 2.2; 2.6; 2.9];
y5meas = [1.75; 2.8; 5.8; 9.3; 12.0; 17.0; 21.0; 25.7];
tmeas  = [1230; 3060; 4920; 7800; 10680; 15030; 22620; 36420];

%Build OPTI Object
theta0 = [0;0;0;0;0]; %inital parameter guess
opts = optiset('display','iter','dynamicOpts',optidynset('sensitivity','cs','initialT',0));
Opt = opti('ode',ode,'data',tmeas,[y1meas' y2meas' y3meas' y4meas' y5meas'],'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% http://tomdyn.com/examples/marinePopulation.html
clc
ode = @(t,z,p) [[0; p(9:end)].*[0; z(1:7)] - (p(1:8)+[p(9:end);0]).*z];

ymeas = [20000 17000 10000 15000 12000 9000 7000 3000
    12445 15411 13040 13338 13484 8426 6615 4022
     7705 13074 14623 11976 12453 9272 6891 5020
     4664  8579 12434 12603 11738 9710 6821 5722
     2977  7053 11219 11340 13665 8534 6242 5695
     1769  5054 10065 11232 12112 9600 6647 7034
      943  3907  9473 10334 11115 8826 6842 7348
      581  2624  7421 10297 12427 8747 7199 7684
      355  1744  5369  7748 10057 8698 6542 7410
      223  1272  4713  6869  9564 8766 6810 6961
      137   821  3451  6050  8671 8291 6827 7525
       87   577  2649  5454  8430 7411 6423 8388
       49   337  2058  4115  7435 7627 6268 7189
       32   228  1440  3790  6474 6658 5859 7467
       17   168  1178  3087  6524 5880 5562 7144
       11    99   919  2596  5360 5762 4480 7256
        7    65   647  1873  4556 5058 4944 7538
        4    44   509  1571  4009 4527 4233 6649
        2    27   345  1227  3677 4229 3805 6378
        1    20   231   934  3197 3695 3159 6454
        1    12   198   707  2562 3163 3232 5566];
tmeas  = 0:0.5:10;

z0 = ymeas(1,:); %ic

%Build OPTI Object
theta0 = [zeros(8,1);zeros(7,1)]; %inital parameter guess
opts = optiset('display','iter','dynamicOpts',optidynset('sensitivity','cs'));
Opt = opti('ode',ode,'data',tmeas,ymeas,'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)


%% http://tomdyn.com/examples/methanolToHydrocarbons.html
clc
ode = @(t,z,p) [-(2*p(2)-(p(1)*z(2))./((p(2)+p(5))*z(1)+z(2))+p(3)+p(4)).*z(1)
                (p(1)*z(1).*(p(2)*z(1)-z(2)))./((p(2)+p(5))*z(1)+z(2))+p(3)*z(1)
                (p(1)*z(1).*(z(2)+p(5)*z(1)))./((p(2)+p(5))*z(1)+z(2))+p(4)*z(1)];

y1meas = [0.7085;0.5971;0.5537;0.3684;0.1712;...
    0.1198;0.0747;0.0529;0.0415;0.0261;0.0208;...
    0.0085;0.0053;0.0019;0.0018];
y2meas = [0.1621;0.1855;0.1989;0.2845;0.3491;...
    0.3098;0.3576;0.3347;0.3388;0.3557;0.3483;...
    0.3836;0.3611;0.3609;0.3485];
y3meas = [0.0811;0.0965;0.1198;0.1535;0.2097;...
    0.2628;0.2467;0.2884;0.2757;0.3167;0.2954;...
    0.295;0.2937;0.2831;0.2846];
tmeas = [0.05;0.065;0.08;0.123;0.233;0.273;...
    0.354;0.397;0.418;0.502;0.553;...
    0.681;0.75;0.916;0.937];

z0 = [1;0;0]; %ic
lb = ones(5,1)*sqrt(eps);
ub = 10*ones(5,1);

%Build OPTI Object
theta0 = [1;1;1;1;1]; %inital parameter guess
opts = optiset('display','iter','solver','nl2sol','dynamicOpts',optidynset('initialT',0,'sensitivity','cs'));
Opt = opti('ode',ode,'bounds',lb,ub,'data',tmeas,[y1meas y2meas y3meas],'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)

%% http://tomdyn.com/examples/parameterEstimation.html
clc
ode = @(t,z,p) [p(1)*z(2); 1-2*z(2)-z(1)];

y1meas = [0.264;0.594;0.801;0.959];
tmeas = [1;2;3;5];

z0 = [NaN;NaN]; %ic

lb = -1.5*ones(3,1); lb(1)=1; %not actually interested in p1, but needed for OPTI
ub = 1.5*ones(3,1); ub(1)=1;

%Build OPTI Object
theta0 = [1;0;0]; %inital parameter guess
opts = optiset('display','iter','solver','nl2sol','dynamicOpts',optidynset('stateIndex',1,'initialT',0,'sensitivity','cs'));
Opt = opti('ode',ode,'bounds',lb,ub,'data',tmeas,y1meas,'theta0',theta0,'z0',z0,'options',opts)

[x,f,e,i] = solve(Opt)
plot(Opt)
