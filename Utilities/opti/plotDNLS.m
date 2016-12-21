function [t,z] = plotDNLS(prob,opts,xb,confStats,tspan)
%plotDNLS Plot Parameter Estimation Problem

%   Copyright (C) 2013 Jonathan Currie (IPL)

if(nargin < 5), tspan = []; end
if(nargin < 4), confStats = []; end

%Measurement Plot Color
measC = [0.4 0.4 0.4];
%Measurement Plot Style
measS = 'o';
%Initial Condition plot style
icS = 'sq';

%Plot confidence lines if present
if(~isempty(confStats))
    plot(NaN,NaN); %hack to prevent patch removing outer borders
    plotConfReg(confStats);
    hold on;
end

% Insert initial conditions if also solved
ind = isnan(prob.odez0);
if(any(ind))
    len = sum(ind);
    estZ0 = xb(end-len+1:end);
    prob.odez0(ind) = estZ0;
    xb = xb(1:end-len);
end

% Generate smooth plot 
dopts = optidynset(opts.dynamicOpts);
if(isempty(tspan))
    tspan = [prob.xdata(1) prob.xdata(end)];
end
ode = @(t,z) prob.ode(t,z,xb);
switch(dopts.integrator)
    case 'ode45'
        [t,z] = ode45(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode15s'
        [t,z] = ode15s(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23'
        [t,z] = ode23(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode113'
        [t,z] = ode113(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23t'
        [t,z] = ode23t(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23tb'
        [t,z] = ode23tb(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode23s'
        [t,z] = ode23s(ode,tspan,prob.odez0,dopts.odeOpts);
    case 'ode15i'
        [t,z] = ode15i(ode,tspan,prob.odez0,dopts.odeOpts);
    otherwise
        optiwarn('OPTI:Plot','Unknown integrator chosen! Using ode45.');
        [t,z] = ode45(ode,tspan,prob.odez0,dopts.odeOpts);
end
if(nargout > 1), return; end %just for evaluating above

%Plot smooth sim
plot(t,z);
title('Dynamic Parameter Estimation Solution');
xlabel('time'); ylabel('z');
hold on;

%Plot Solved Initial Conditions
if(any(ind))
    for i = 1:sum(ind)
        plot(tspan(1),estZ0(i),icS,'Color',measC);
    end
end

%If we have xdata_old, then should be in cell format, plot each cell
if(isfield(prob.misc,'xdata_orig') && ~isempty(prob.misc.xdata_orig))
    if(iscell(prob.misc.xdata_orig) && iscell(prob.misc.ydata_orig))
        for i = 1:length(prob.misc.xdata_orig)
            plot(prob.misc.xdata_orig{i},prob.misc.ydata_orig{i},measS,'Color',measC);
        end 
    elseif(~iscell(prob.misc.xdata_orig) && ~iscell(prob.misc.ydata_orig))
        for i = 1:size(prob.misc.ydata_orig,2)
            plot(prob.misc.xdata_orig,prob.misc.ydata_orig(:,i),measS,'Color',measC);  
        end
    else
        error('Expected misc.xdata_orig and misc.ydata_orig to be cell arrays');
    end    
else
    %Otherwise reshape ydata (always a column in OPTI) and plot
    if(~isempty(dopts.stateIndex))
        nstates = length(prob.odez0(dopts.stateIndex));
    else
        nstates = length(prob.odez0);
    end
    ydata = reshape(prob.ydata,length(prob.xdata),nstates);
    %Plot
    plot(prob.xdata,ydata,measS,'Color',measC);
end

hold off;