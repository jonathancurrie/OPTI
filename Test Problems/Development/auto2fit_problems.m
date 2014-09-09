%% Auto2Fit Test Problems
%http://www.7d-soft.com/en/AF_Examples.htm#5:%20Chart
clear all
clc
load auto2fitData

%% Problem 1.1 - Mount Curve
clc
%Fitting Function
fun = @(m,xdata) m(1)/m(2) * exp(-0.5 * ((xdata-m(3))./m(2)).^2);
%Data
xdata = data11(:,2);
ydata = data11(:,3);

%Problem Setup
opts = optiset('display','iter','solver','nomad');
Opt = opti('fun',fun,'data',xdata,ydata,'x0',[1;3;460],'bounds',[0;0;100],[5;10;1000],'options',opts)

%Solve & Plot
[x,f,e,i] = solve(Opt)
plot(Opt)


%% Problem 1.1 - Mount Curve [multi-solve]
clc
%Fitting Function
fun = @(m,xdata) m(1)/m(2) * exp(-0.5 * ((xdata-m(3))./m(2)).^2);
%Data
xdata = data11(:,2);
ydata = data11(:,3);

%Problem Setup
opts = optiset('display','final','solver','nl2sol');
Opt = opti('fun',fun,'data',xdata,ydata,'bounds',[0;0;100],[5;10;1000],'options',opts)

%Solve & Plot
[x,f,e,i] = multisolve(Opt,[],12)
plot(Opt)


%% Problem 1.2 - Cos Function
clc
%Fitting Function
fun = @(m,xdata) (m(1)./xdata(:,1) - cos(m(2).*xdata(:,2))) .* xdata(:,3)./xdata(:,1) .* m(3);
%Data
xdata = data12(:,2:4);
ydata = data12(:,end);

%Problem Setup
opts = optiset('display','iter','solver','nomad');
Opt = opti('fun',fun,'data',xdata,ydata,'x0',[1;40;1],'bounds',[0;0;0],[5;100;5],'options',opts)

%Solve & Plot
[x,f,e,i] = solve(Opt)
plot(Opt)


%% Problem 1.2 - Cos Function [multi-solve]
clc
%Fitting Function
fun = @(m,xdata) (m(1)./xdata(:,1) - cos(m(2).*xdata(:,2))) .* xdata(:,3)./xdata(:,1) .* m(3);
%Data
xdata = data12(:,2:4);
ydata = data12(:,end);

%Problem Setup
opts = optiset('display','final','solver','nl2sol');
Opt = opti('fun',fun,'data',xdata,ydata,'x0',[1;40;1],'bounds',[0;0;0],[5;100;5],'options',opts)

%Solve & Plot
[x,f,e,i] = multisolve(Opt)
plot(Opt)


%% Problem 1.2 - Cos Function [multi-solve, 2 variable version]
clc
%Fitting Function
fun = @(m,xdata) (2.4923./xdata(:,1) - cos(m(1)*xdata(:,2))) .* xdata(:,3)./xdata(:,1) .* m(2);
%Data
xdata = data12(:,2:4);
ydata = data12(:,end);

%Problem Setup
opts = optiset('display','final','solver','auto');
Opt = opti('fun',fun,'data',xdata,ydata,'bounds',[0;0],[100;5],'options',opts)

%Solve & Plot
[x,f,e,i] = multisolve(Opt,[],[5 3])

figure(1);
plot(Opt)
figure(2);
multiplot(Opt,1,75)


%% Problem 1.3 - SinCos Function (data appears wrong)
% clc
% %Fitting Function
% fun = @(a,xdata) (a(1) + xdata(:,1)./a(2) + cos(a(3) * xdata(:,2)./xdata(:,3)))./(a(4) .* sin(xdata(:,1) + xdata(:,2) + xdata(:,3)));
% %Data
% xdata = data13(:,2:4);
% ydata = data13(:,end);
% 
% %Problem Setup
% opts = optiset('display','iter','solver','nomad');
% Opt = opti('fun',fun,'data',xdata,ydata,'x0',[1 1 10 1],'options',opts)
% 
% %Solve & Plot
% [x,f,e,i] = solve(Opt)
% plot(Opt)

