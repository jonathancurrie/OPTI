function genConfidenceData
%% Generate unit test data for confidence limits from MATLAB's statistics toolbox
% J.Currie October 2017

% Himmelblau Example
% Function r = f(x,theta) (x = pressure, theta unknowns)
Rxn_rate = @(theta,p) theta(1)*p./(1+theta(2)*p); % r = f(x,theta) (x = pressure, theta unknowns)
% Fitting Data
p   = [20,30,35,40,50,55,60]'; % Pressure, P
r   = [0.068,0.0858,0.0939,0.0999,0.1130,0.1162,0.1190]'; % Reaction rate, r
x0  = [5e-3 2e-2];
wts = ones(size(p));

himmel  = fitModel(Rxn_rate, p, r, x0, wts);


% Himmeblau with Weights
wts(end) = 10;
himmelw  = fitModel(Rxn_rate, p, r, x0, wts);


% SAS Example
% Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
% Fitting Data
n   = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r   = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r
x0  = [150 0.01];
wts = ones(size(n));

sas  = fitModel(eKin, n, r, x0, wts);

% SAS Example with weights
wts(end) = 10;
sasw  = fitModel(eKin, n, r, x0, wts);


% Save to mat
clear Rxn_rate p r x0 wts eKin n
save Utilities/UnitTests/UnitTestData/confidenceData




% Save parameters from nonlinear model to generic struct
function data = fitModel(fun, xdata, ydata, x0, wts)

mdl = NonLinearModel.fit(xdata,ydata,fun,x0,'weights',wts);

cis = coefCI(mdl); 
[~,ys] = predict(mdl,xdata,'Simultaneous','true');

data.xdata      = xdata;
data.ydata      = ydata;
data.x0         = x0;
data.wts        = wts;
data.Rsquare    = mdl.Rsquared.Ordinary;
data.AdjRsquare = mdl.Rsquared.Adjusted;
data.RMSE       = mdl.RMSE;
data.DFE        = mdl.DFE;
if (isa(mdl.Coefficients,'table'))
    data.param      = table2array(mdl.Coefficients);
else
    data.param      = mdl.Coefficients;
end
data.confInt    = cis(:,2) - mdl.Coefficients.Estimate;
data.confBnds   = ys;

