%% OPTIFIT Testing
clc
clear all

%% EXP1
clc
%True Function
fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata);
theta = [400;-0.1];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'exp1')
plot(ofit)

%% EXP2 +
clc
%Function
fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3);
theta = [2;0.1;5];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'exp2')
plot(ofit)

%% EXP2 + w NEG
clc
%Function
fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3);
theta = [2;0.1;-3000];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'exp2')
plot(ofit)

%% EXP2 -
clc
%Function
fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3);
theta = [2;-0.1;5];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'exp2')
plot(ofit)

%% EXP3
clc
%Function
fun = @(theta,xdata) theta(1)*exp(theta(2)*xdata) + theta(3)*exp(theta(4)*xdata);
theta = [2;-0.1;50;-0.05];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'exp3')
plot(ofit)

%% POWER1
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^theta(2);
theta = [3;0.5];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'power1')
plot(ofit)

%% POWER2 +
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^theta(2) + theta(3);
theta = [3;0.5;-100];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'power2')
plot(ofit)

%% POWER2 -
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^theta(2) + theta(3);
theta = [3;-0.5;-100];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'power2')
plot(ofit)


%% NLS2
clc
%True Function
xdata = (1:40)';
fun = @(x,xdata) x(1)*exp(-x(2)*xdata) + x(3); %actual
%Fitting Data
ydata=[5.8728, 5.4948, 5.0081, 4.5929, 4.3574, 4.1198, 3.6843, 3.3642, 2.9742, 3.0237, 2.7002, 2.8781,...
       2.5144, 2.4432, 2.2894, 2.0938, 1.9265, 2.1271, 1.8387, 1.7791, 1.6686, 1.6232, 1.571, 1.6057,...
       1.3825, 1.5087, 1.3624, 1.4206, 1.2097, 1.3129, 1.131, 1.306, 1.2008, 1.3469, 1.1837, 1.2102,...
       0.96518, 1.2129, 1.2003, 1.0743];

ofit = optifit(xdata,ydata,'auto')
plot(ofit)


%% POLY TEST
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^theta(2) + theta(3);
theta = [3;0.5;-100];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'poly1')
plot(ofit)

%% RAT01
clc
%Function
fun = @(theta,xdata) theta(1) ./ (xdata + theta(2));
theta = [3;0.5];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat01')
plot(ofit)

%% RAT02
clc
%Function
fun = @(theta,xdata) theta(1) ./ (xdata.^2 + theta(2)*xdata + theta(3));
theta = [3;0.5;7];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat02')
plot(ofit)

%% RAT03
clc
%Function
fun = @(theta,xdata) theta(1) ./ (xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4));
theta = [3;0.5;7;1];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat03')
plot(ofit)

%% RAT11
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata + theta(3));
theta = [3;0.5;2];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat11')
plot(ofit)

%% RAT12
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata.^2 + theta(3)*xdata + theta(4));
theta = [3;0.5;2;5];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat12')
plot(ofit)

%% RAT13
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata + theta(2)) ./ (xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5));
theta = [3;0.5;2;5;0.1];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3]';
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat13')
plot(ofit)

%% RAT21
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata + theta(4));
theta = [0.1;0.3;5.8;0.06];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat21')
plot(ofit)

%% RAT22
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata.^2 + theta(4)*xdata + theta(5));
theta = [0.1;0.3;5.8;0.06;1];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat22')
plot(ofit)

%% RAT23
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata.^3 + theta(4)*xdata.^2 + theta(5)*xdata + theta(6));
theta = [0.1;0.3;5.8;0.06;1;-0.1];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'rat23')
plot(ofit)

%% SIN1
clc
fcn = @(p,t) p(1)*sin(p(2)*t+p(3)) + 5;
px = [220 2*pi*50 0.5]'; 
t = linspace(0,1/20,80)'; 
y = fcn(px,t) + 50*randn(size(t)); 

ofit = optifit(t,y,'sin1')
plot(ofit)

%% SIN2
clc
clear
%Sum of Sine Problem
fcn = @(p,t) p(1)*sin(p(2)*t+p(3)) + p(4)*sin(p(5)*t+p(6)) + p(7);
px = [220 2*pi*53.5 0.5, 120 2*pi*45 0.2 100]'; 
t = linspace(0,1/10,100)';
% t = sort(t(randi(length(t),ceil(0.75*length(t)),1))); %random samples (makes it very tricky)
y = fcn(px,t) + 10*randn(size(t)); 

ofit = optifit(t,y,'sin2')
plot(ofit)


%% SIN3
clc
clear
%Sum of Sine Problem
fcn = @(p,t) p(1)*sin(p(2)*t+p(3)) + p(4)*sin(p(5)*t+p(6)) + p(7)*sin(p(8)*t+p(9)) + p(10);
px = [220 2*pi*53.5 0.5, 120 2*pi*45 0.2, 60 2*pi*25 0.8 30]'; 
t = linspace(0,1/10,100)';
% t = sort(t(randi(length(t),ceil(0.75*length(t)),1))); %random samples (makes it very tricky)
y = fcn(px,t) + 10*randn(size(t)); 

ofit = optifit(t,y,'sin3')
plot(ofit)


%% RAT21 AUTO
clc
%Function
fun = @(theta,xdata) (theta(1)*xdata.^2 + theta(2)*xdata + theta(3)) ./ (xdata + theta(4));
theta = [0.1;0.3;5.8;0.06];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'auto')
plot(ofit)

%% POLY3 AUTO
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata + theta(4);
theta = [0.05;-3;-0.6;0.06];
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'auto')
plot(ofit)

%% POLY4 AUTO
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^4 + theta(2)*xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5);
theta = [0.05;-3;-6;0.06;1]*0.0001;
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);

ofit = optifit(xdata,ydata,'auto')
plot(ofit)

%% SIN2 AUTO
clc
fcn = @(p,t) p(1)*sin(p(2)*t+p(3)) + p(4)*sin(p(5)*t+p(6)) + p(7);
px = [220 2*pi*53.5 0.5, 120 2*pi*45 0.2 100]'; 
t = linspace(0,1/10,100)';
y = fcn(px,t);

ofit = optifit(t,y,'auto')
plot(ofit)


%% 3D POLY11
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata + theta(2)*ydata + theta(3);
theta = [3;0.5;3];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly11')
plot(ofit)

%% 3D POLY12
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*ydata.^2 + theta(2)*xdata.*ydata + theta(3)*xdata + theta(4)*ydata + theta(5);
theta = [3;0.5;3;2;1];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly12')
plot(ofit)

%% 3D POLY13
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*ydata.^3 + theta(2)*ydata.^2 + theta(3)*xdata.*ydata.^2 + theta(4)*xdata.*ydata + theta(5)*xdata + theta(6)*ydata + theta(7);
theta = [3;0.5;3;2;1;0.5;2];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly13')
plot(ofit)

%% 3D POLY21
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^2 + theta(2)*xdata.*ydata + theta(3)*xdata + theta(4)*ydata + theta(5);
theta = [3;0.5;3;2;1];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly21')
plot(ofit)

%% 3D POLY22
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^2 + theta(2)*ydata.^2 + theta(3)*xdata.*ydata + theta(4)*xdata + theta(5)*ydata + theta(6);
theta = [3;0.5;3;2;1;7];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly22')
plot(ofit)


%% 3D POLY23
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*ydata.^3 + theta(2)*xdata.^2 + theta(3)*ydata.^2 + theta(4)*xdata.^2.*ydata + theta(5)*xdata.*ydata.^2 + theta(6)*xdata.*ydata + theta(7)*xdata + theta(8)*ydata + theta(9);
theta = [3;0.5;3;2;1;0.5;2;1;4];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly23')
plot(ofit)

%% 3D POLY31
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*xdata.^2.*ydata + theta(4)*xdata.*ydata + theta(5)*xdata + theta(6)*ydata + theta(7);
theta = [3;0.5;3;2;1;-1;8];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly31')
plot(ofit)

%% 3D POLY32
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*xdata.^2 + theta(3)*ydata.^2 + theta(4)*xdata.^2.*ydata + theta(5)*xdata.*ydata.^2 + theta(6)*xdata.*ydata + theta(7)*xdata + theta(8)*ydata + theta(9);
theta = [3;0.5;3;2;1;-1;8;-5;10];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly32')
plot(ofit)

%% 3D POLY33
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*ydata.^3 + theta(3)*xdata.^2 + theta(4)*ydata.^2 + theta(5)*xdata.^2.*ydata + theta(6)*xdata.*ydata.^2 + theta(7)*xdata.*ydata + theta(8)*xdata + theta(9)*ydata + theta(10);
theta = [3;0.5;3;2;1;-1;8;-5;10;1];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},'poly33')
plot(ofit)

%% 2D General Fun
clc
%Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
%Fitting Data
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r

ofit = optifit(n,r,eKin,[0 0])
plot(ofit)

%% 3D General Fun
clc
%Function
fun = @(theta,xdata,ydata) theta(1)*xdata.^3 + theta(2)*ydata.^3 + theta(3)*xdata.^2 + theta(4)*ydata.^2 + theta(5)*xdata.^2.*ydata + theta(6)*xdata.*ydata.^2 + theta(7)*xdata.*ydata + theta(8)*xdata + theta(9)*ydata + theta(10);
theta = [3;0.5;3;2;1;-1;8;-5;10;1];
%Fitting Data
xdata = (1:10)';
ydata = (1:10)';
[X,Y] = meshgrid(xdata,ydata);
Z = fun(theta,X,Y);

ofit = optifit(X,{Y Z},fun,zeros(size(theta)))
plot(ofit)

%% SAS Example
clc
%Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
%Fitting Data
n = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r

ofit = optifit(n,r)
plot(ofit)

%% Stats Toolbox Example
clc
if (exist('NonlinearModel.m','file'))
    load carbig
    X = Horsepower;
    Y = Weight;
    Z = MPG;
    modelfun = @(b,xdata,ydata)b(1) + b(2)*xdata.^b(3) + b(4)*ydata.^b(5);
    beta0 = [-50 500 -1 500 -1];

    ofit = optifit(X,{Y Z},modelfun,0*beta0)
    plot(ofit)
end

%% Himmelblau
clc
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r

ofit = optifit(p,r)
plot(ofit)

%% Himmelblau NL
clc
p = [20,30,35,40,50,55,60]'; % Pressure, P
r = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r

ofit = optifit(p,r)
plot(ofit)
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)

%% POLY4 w WTS
clc
%Function
fun = @(theta,xdata) theta(1)*xdata.^4 + theta(2)*xdata.^3 + theta(3)*xdata.^2 + theta(4)*xdata + theta(5);
theta = [0.05;-3;-6;0.06;1]*0.0001;
%Fitting Data
xdata = [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = fun(theta,xdata);% + randn(size(xdata));
wts = ones(size(ydata));
wts(3) = 100;

ofit = optifit(xdata,ydata,'poly3',[],[],[],wts)
plot(ofit)