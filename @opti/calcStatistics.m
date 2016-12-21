function stats = calcStatistics(optObj,limit)
%Calculate confidence and variance statistics of an OPTI curve fitting problem
%
%   Called By opti fitStats

%   Copyright (C) 2011-2014 Jonathan Currie (IPL)

% This function uses ideas from:
% - David M. Himmelblau. Process Analysis by Statistical Methods. 
%   John Wiley & Sons, 1970
% - David I. Wilson. Advanced Control Using MATLAB, AUT University, 2014.
% - NonLinearModel.m from the Statistics Toolbox
% - confint.m from the Curve Fitting Toolbox
% - predint.m from the Curve Fitting Toolbox

% Value explanations
% http://www.ats.ucla.edu/stat/sas/output/reg.htm

if(nargin < 2 || isempty(limit)), limit = 0.95; end

if(limit <= 0 || limit >= 1)
    error('Confidence limit must be 0 < limit < 1');
end

%Get Problem
prob = optObj.prob;
x = optObj.sol;
nparam = length(x);
ndata = length(prob.ydata); %check whether we need original data?
wlevel = optiWarnLevel(optObj.opts.warnings);
%Get Solution SSE
sse = optObj.obj;
dfe = ndata-nparam;
rmse = sqrt(sse/dfe);
if(wlevel>1)
    optiwarn('optistats:dfe','Your problem has #data <= #parameters, confidence limits and most statistics cannot be estimated');
end
%Build Return Structure
stats = struct('SSE',sse,'Rsquare',NaN,'AdjRsquare',NaN,'RMSE',NaN,'DFE',ndata-nparam,'ConfInt',NaN,'ConfBnds',struct('xdata',NaN,'bnds',NaN),'Cov',NaN,'BndIdx',NaN,...
               'Param',struct('StdError',NaN,'tStat',NaN,'pValues',NaN),'Model',struct('FStat',NaN,'pValue',NaN),'Conf',limit,'SolverStatus','');

nidx = strfind(prob.Name,'optifit:');           
if(~isempty(nidx))           
    stats.ModelStructure = prob.Name(9:end);
end

%Don't be overly enthusiastic
if(~isempty(strfind(optObj.info.Status,'Optimal')))
    stats.SolverStatus = 'Converged';
else
    stats.SolverStatus = optObj.info.Status;
end
           
%Ensure we have a NLS or DNLS
if(~any(strcmpi(prob.type,{'NLS','DNLS'})))
    error('This function can only be used on curve fitting (NLS,DNLS) problems');
end

%Get ydata and weights
ydata = prob.ydata;
weights = prob.weighting;

%Get Partial Derivative Matrix (data/parameters) [original xdata]
if(~isempty(prob.xdata))
    try
        X = prob.misc.fitGrad(x,prob.xdata);        
        gmode = 2;
    catch
        X = prob.misc.fitGrad(x);        
        gmode = 1;
    end
else
    try
        X = prob.misc.fitGrad(x);
    catch
        error('Unknown gradient field!');
    end
    gmode = 1;
end
%If we have weights, weight the Jacobian
if(~isempty(weights))
    X = bsxfun(@times, weights, X);
end
%Check for active bounds, if so, remove cols of X
bidx = false(size(X,2),1);
if(~isempty(prob.lb) || ~isempty(prob.ub))
    if(any(~isinf(prob.lb))) 
        bidx = bidx | abs(x-prob.lb) < sqrt(eps);
    end
    if(any(~isinf(prob.ub))) 
        bidx = bidx | abs(x-prob.ub) < sqrt(eps);
    end
end
if(any(bidx))
    if(wlevel > 1)
        optiwarn('opti:conf','One or more bounds are active at the solution. The confidence interval and parameter statistics cannot be calculated for these variables.');
    end
    X = X(:,~bidx);
    %Correct dfe, rmse
    nparam = nparam - sum(bidx);
    dfe = ndata-nparam;
    rmse = sqrt(sse/dfe);
end

try
    %Solve inverse of X'*X, but only interested in diagonal
    [~,R] = qr(X,0);
    if(size(R,1)~=size(R,2))
        throw(MException('opti:conf','QR R not square'));
    end
    s = warning('off','MATLAB:singularMatrix');
    s1 = warning('off','MATLAB:nearlySingularMatrix');
    Rinv = R \ eye(length(R));
    warning(s1);
    warning(s);
catch
    error('The problem is poorly scaled and resulted in a (near) singular matrix, limits cannot be calculated');
end

%Solve Covariance Matrix
s = warning('off','MATLAB:singularMatrix');
s1 = warning('off','MATLAB:nearlySingularMatrix');
try
    L = chol(X'*X,'lower');
    stats.Cov = L'\(L\eye(length(L)))*rmse^2;
catch
    if(wlevel)
        optiwarn('opti:conf','The X''*X matrix has been found to be singular, and a general inverse will be used to (attempt to) solve the system covariance.\nResults are now suspicious at best...');
    end    
    stats.Cov = inv(X'*X)*rmse^2; %poor method - any better ideas?    
end
warning(s1);
warning(s);
stats.BndIdx = bidx;

%Generate ypred for ANOVA stuff
if(nargin(prob.misc.fitFun)==2)
    sumypred = sum(prob.misc.fitFun(x,prob.xdata).^2);
else
    sumypred = sum(prob.misc.fitFun(x).^2);
end

%Solve Confidence Interval
v = sum(Rinv.^2,2) * (sse / dfe);
alpha = (1-limit)/2;
ConfInt = (-rmathlib('qt',alpha,dfe)* sqrt(v'))';
stats.ConfInt = NaN(length(bidx),1);
stats.ConfInt(~bidx) = ConfInt;

%Solve Confidence Bounds for Plotting (note we must have system covariance for this to work!)
if(~isempty(stats.Cov))
    %If we have an ODE and simple time step (i.e. not a cell array or repeated measurements, try for a smooth curve)
    if(isfield(prob.misc,'xdata_orig') && ~isempty(prob.misc.xdata_orig) && ~iscell(prob.misc.xdata_orig))
        havRep = length(unique(prob.misc.xdata_orig)) ~= length(prob.misc.xdata_orig);
    else
        havRep = false;
    end
    if(~isempty(prob.ode) && ~iscell(prob.misc.xdata_orig) && length(unique(prob.xdata))==length(prob.xdata) && ~havRep)
        try
            %Lazy way to convert problem again
            xd = linspace(min(prob.xdata),max(prob.xdata),max(10*length(prob.xdata),1e2))';
            prob.ydata = ones(length(prob.ydata)/length(prob.xdata)*length(xd),1); prob.xdata = xd; 
            prob = DNLS2NLS(prob,optObj.opts);
            Xb = prob.misc.fitGrad(x);
            ypred = prob.misc.fitFun(x,xd);
        catch %no luck
            xd = prob.xdata;
            Xb = X;
            if(nargin(prob.misc.fitFun)==2)
                ypred = prob.misc.fitFun(x,xd);    
            else
                ypred = prob.misc.fitFun(x);
            end
        end
    else %algebraic system or complicated ode problem
        if(gmode==2)
            %try evaluate at many intermediate points for a smoother curve
            try
                %Assume if xdata is a matrix, or not sorted, don't try smooth
                if(size(prob.xdata,1) > 1 && size(prob.xdata,2) > 1 || (any(sort(prob.xdata) ~= prob.xdata) && any(sort(prob.xdata,'descend') ~= prob.xdata)))
                    error('skip');
                end        
                xd = linspace(min(prob.xdata),max(prob.xdata),max(10*length(prob.xdata),1e2))';
                Xb = prob.misc.fitGrad(x,xd);
                ypred = prob.misc.fitFun(x,xd);
            catch
                xd = prob.xdata;
                Xb = X;
                ypred = prob.misc.fitFun(x,prob.xdata);
            end
        else %no luck
            xd = prob.xdata;
            Xb = X;
            if(nargin(prob.misc.fitFun)==2)
                %if function has two input arguments, try a numerical gradient
                try
                    if(~isempty(strfind(char(prob.misc.fitFun),'odeEstim'))), error('no luck'); end %integrator most likely has set time points it is expecting
                    xd = linspace(min(prob.xdata),max(prob.xdata),max(10*length(prob.xdata),1e2))';
                    Xb = mklJac(@(x) prob.misc.fitFun(x,xd),x);
                    ypred = prob.misc.fitFun(x,xd);
                catch
                    xd = prob.xdata;
                    Xb = X;
                    ypred = prob.misc.fitFun(x,xd);
                end
            else
                ypred = prob.misc.fitFun(x);
            end
        end
    end

    %Ensure ypred is a column
    if(size(ypred,2) > 1), ypred = ypred'; end
    %Solve Prediction Bounds (Simultaneous Functional)    
    crit = sqrt(length(x) * rmathlib('qf',limit, length(x), dfe));
    %If we have active bounds, best to calculate Jacobian at each point we want to plot, then drop columns, rather than use Covariance below
    if(any(bidx))
        try
            J = prob.misc.fitGrad(x,xd); 
            E = J(:,~bidx)*Rinv;
            delta = crit * sqrt(sum(E.*E,2)) * sqrt(sse/dfe);           
        catch
            delta = crit * sqrt(sum((Xb*stats.Cov) .* Xb,2));
        end
    else
        delta = crit * sqrt(sum((Xb*stats.Cov) .* Xb,2));
    end
    if(isempty(xd))
        stats.ConfBnds.xdata = (1:length(ypred))';
    else
        stats.ConfBnds.xdata = xd;
    end
    stats.ConfBnds.bnds = [ypred-delta ypred+delta];
    

    %Reshape based on number of curves
    if(~isempty(prob.odez0))
        if(~isempty(optObj.opts.dynamicOpts) && isfield(optObj.opts.dynamicOpts,'stateIndex') && ~isempty(optObj.opts.dynamicOpts.stateIndex))
            nc = length(optObj.opts.dynamicOpts.stateIndex); 
        else
            nc = length(prob.odez0);
        end    
        nd = size(ypred,1)/nc;
        if(floor(nd)==nd) %have to work out xdata associated with which point... anyone interested?
            cb = zeros(nd,nc*2);
            idx = 1;
            for i = 1:2:nc*2
                cb(:,i:i+1) = stats.ConfBnds.bnds(idx:idx+nd-1,:);
                idx = idx + nd;
            end
            stats.ConfBnds.bnds = cb;
        else
            stats.ConfBnds.bnds = NaN;
            stats.ConfBnds.xdata = NaN;
            if(wlevel)
                optiwarn('opti:confbnds','Cannot currently determine confidence prediction bounds for this type of problem, sorry!'); 
            end
        end
    end
    %If we had repeated points, try index out...
    if(havRep)
        %Check for repeated measurements
        try
            [~,ia] = unique(prob.misc.xdata_orig,'legacy');
        catch
            [~,ia] = unique(prob.misc.xdata_orig);
        end
        stats.ConfBnds.bnds = stats.ConfBnds.bnds(ia,:);
    end
    
    %Make sure in order for plotting
    [stats.ConfBnds.xdata,sbidx] = sort(stats.ConfBnds.xdata);
    stats.ConfBnds.bnds = stats.ConfBnds.bnds(sbidx,:);
end

%Solve Standard Statistics
[SST,stats.Rsquare,stats.AdjRsquare,stats.RMSE,SST0] = fitstats(ydata,weights,sse,dfe,ndata);
smodel = sumypred-sse;
stats.DF = [nparam dfe ndata]; %[model error uncorrected]
stats.SOS = [smodel sse sumypred];
stats.MS = [smodel/nparam sse/(ndata-nparam)]; 

%Parameter statistics require covariance too...
if(~isempty(stats.Cov))
    %Parameter Statistics
    SE = NaN(nparam,1); T = SE; P = SE;
    SE(~bidx) = sqrt(diag(stats.Cov));
    T(~bidx) = x(~bidx) ./ SE(~bidx);
    P(~bidx) = 2*rmathlib('pt',-abs(T(~bidx)),dfe);
    stats.Param.StdError = SE;
    stats.Param.tStat = T;
    stats.Param.pValues = P;
end

%Model F-Test & P Value
%Check for intercept (constant column in X)
Xmin = min(X,[],1); Xmax = max(X,[],1);
if(nparam > 1 && any(abs(Xmax-Xmin) <= sqrt(eps)))
    %Intercept Model    
    DFR = nparam - 1;
    DFE = ndata - 1 - DFR;
    SSR = max(SST - sse,0);
    stats.Model.FStat = (SSR/DFR) / (sse / DFE);
    stats.Model.pValue = rmathlib('pf',1/stats.Model.FStat,DFE,DFR);
    stats.Model.Int = true;
else
    %Zero Model
    DFR = nparam;
    DFE = ndata - nparam;
    SSR = max(SST0 - sse,0);
    stats.Model.FStat = (SSR/DFR) / (sse / DFE);
	stats.Model.pValue = rmathlib('pf',1/stats.Model.FStat,DFE,DFR);
    stats.Model.Int = false;
end


function [sst,rsquare,adjrsquare,rmse,sst0] = fitstats(ydata,weights,sse,dfe,ndata)
%Calculate standard statistical measures

%Spread about mean
if(isempty(weights))
    ybar = mean(ydata);
    sst = sum((ydata - ybar).^2);
    sst0 = sum(ydata.^2);
else
    ybar = sum(ydata.*weights.^2)/sum(weights.^2);
    sst = sum(weights.^2.*(ydata - ybar).^2);
    sst0 = sum(weights.^2.*ydata.^2);
end
%R^2
if(sst~=0)
    rsquare = 1 - sse/sst;
    adjrsquare = 1 - (1-rsquare)*(ndata-1)/dfe; %1 - (sse/dfe)/(sst/ndata);
else
    rsquare = NaN;
    adjrsquare = NaN;
end
%RMSE
if(dfe > 0)
    rmse = sqrt(sse/dfe);
else
    rmse = NaN;
end

