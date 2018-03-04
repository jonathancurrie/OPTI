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

% Himmelblau with Weights
wts(end) = 10;
himmelw  = fitModel(Rxn_rate, p, r, x0, wts);

% Himmelblau with lb active
lb = [0;0.03];
ub = [7e-3,1];
wts = ones(size(p));
ft = fittype( 'a*x./(1+b*x)', 'independent', 'x', 'dependent', 'y' );
himmellb = fitModelBounded(ft, p, r, x0, lb, ub, wts);

% Himmelblau with ub active
lb = [0;0];
ub = [5e-3,1];
himmelub = fitModelBounded(ft, p, r, x0, lb, ub, wts);

% Himmelblau with both bounds
lb = [0;0.03];
ub = [5e-3,1];
himmelbnd = fitModelBounded(ft, p, r, x0, lb, ub, wts);

% Himmelblau with bounds and weights
lb = [0;0];
ub = [4e-3,1];
wts(end) = 10;
himmelbndw = fitModelBounded(ft, p, r, x0, lb, ub, wts);


% SAS Example
% Function
eKin = @(theta,n) theta(1)*n./(theta(2) + n);
% Fitting Data
n   = [0.26,0.3,0.48,0.5,0.54,0.64,0.82,1.14,1.28,1.38,1.8,2.3,2.44,2.48]'; % Amount of Substrate
r   = [124.7,126.9,135.9,137.6,139.6,141.1,142.8,147.6,149.8,149.4,153.9,152.5,154.5,154.7]'; % Reaction rate, r
x0  = [150 0.01];
wts = ones(size(n));
sas = fitModel(eKin, n, r, x0, wts);

% SAS Example with weights
wts(end) = 10;
sasw  = fitModel(eKin, n, r, x0, wts);


% Save to mat
clear Rxn_rate p r x0 wts eKin n lb ub
save Utilities/UnitTests/UnitTestData/confidenceData




% Save parameters from nonlinear model to generic struct
function data = fitModel(fun, xdata, ydata, x0, wts)

mdl = NonLinearModel.fit(xdata,ydata,fun,x0,'weights',wts);

cis = coefCI(mdl); 
[~,ys] = predict(mdl,xdata,'Simultaneous','true','Weights',wts);

data.xdata      = xdata;
data.ydata      = ydata;
data.x0         = x0;
data.wts        = wts;
data.lb         = [];
data.ub         = [];
data.Rsquare    = mdl.Rsquared.Ordinary;
data.AdjRsquare = mdl.Rsquared.Adjusted;
data.RMSE       = mdl.RMSE;
data.DFE        = mdl.DFE;
if (isa(mdl.Coefficients,'table'))
    data.param  = table2array(mdl.Coefficients);
    data.sol    = mdl.Coefficients.Estimate;
else
    data.param  = mdl.Coefficients;
    data.sol    = mdl.Coefficients(:,1);
end
data.confInt    = cis(:,2) - mdl.Coefficients.Estimate;
data.confBnds   = ys;
data.SSE        = mdl.SSE;
data.cov        = mdl.CoefficientCovariance;


% Save paramters from curve fitting toolbox model to generic struct
function data = fitModelBounded(fun, xdata, ydata, x0, lb, ub, wts)

[xData, yData, weights] = prepareCurveData( xdata, ydata, wts );
opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display    = 'Off';
opts.Lower      = lb;
opts.StartPoint = x0;
opts.TolFun     = 1e-09;
opts.TolX       = 1e-09;
opts.Upper      = ub;
opts.Weights    = weights;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, fun, opts );

data.xdata      = xdata;
data.ydata      = ydata;
data.x0         = x0;
data.wts        = wts;
data.lb         = lb;
data.ub         = ub;
data.Rsquare    = gof.rsquare;
data.AdjRsquare = gof.adjrsquare;
data.RMSE       = gof.rmse;
data.DFE        = gof.dfe;
data.sol        = coeffvalues(fitresult);
data.sol        = data.sol(:);
ci              = confint(fitresult)';
data.confInt    = ci(:,2) - data.sol;
if (~all(isnan(data.confInt)))
    data.confBnds = predint(fitresult,xdata,0.95,'Functional','on');
else
    data.confBnds = NaN;
end
data.SSE        = gof.sse;




